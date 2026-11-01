# GimmeMotifs：基序扫描、聚类、背景集合与从头发现

这一组 40 个函数是 `dnakit.sequence` 中对 [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs)（转录因子基序分析工具集）的适配器。它们用位置频率矩阵（PFM）扫描序列，对基序做匹配、聚类和绘图，构建背景序列集合（随机、基因组和 GC 匹配），合并峰 summit，统计峰内 reads 数，比较序列集合之间的基序频率，写出 ROC 与富集报告，运行从头（de novo）基序预测工具，并根据同源群构建基序到因子的对应表。部分函数还会用到 [pybedtools](https://github.com/daler/pybedtools) 与 [bedtools](https://github.com/arq5x/bedtools2) 程序、[ete3](https://github.com/etetoolkit/ete)、[gffread](https://github.com/gpertea/gffread) 和 [OrthoFinder](https://github.com/davidemms/OrthoFinder)。所有函数都处理本地文件，不会下载基因组或数据库；基因组始终以 FASTA 路径给出，而不是 genomepy 名称。本组中唯一访问网络的函数是 `gimmemotifs_jaspar_motif_info()`。

GimmeMotifs 使用 MIT 许可（已安装包中的许可文件是 MIT 许可）。这些函数从 `dnakit.sequence` 导入，例如 `from dnakit.sequence import gimmemotifs_pfm_scan`。依赖都是可选的，只有调用对应函数时才会导入：用 `pip install "dnakit[bio-gimmemotifs]"` 安装（GimmeMotifs 0.18 与 ete3）。外部程序（`gimme`、`bedtools`、`gffread`、`orthofinder` 以及各从头预测工具）在函数运行时才查找，找不到就抛出 `BackendUnavailableError`。除 `gimmemotifs_jaspar_motif_info()` 外，每个函数都返回一个字典，含 `status`、`operation`、`backend`、`version` 和 `outputs`（输出名对应路径字符串），外加各节中列出的字段；`gimmemotifs_jaspar_motif_info()` 只返回 `motif_id` 和 `info`。输出文件和目录必须是新的：已存在的路径或目录会被 `ConfigurationError` 拒绝。唯一的例外是 `gimmemotifs_maelstrom_html_report()`，它把报告写进已存在的运行目录。

## 示例数据

大多数示例使用 `motifs.pfm` 中的三个小基序、`seqs.fa` 中的六条短序列（其中部分序列含 `TGACAC` 或 `GCATAC` 位点），以及 `db.pfm` 中的一个数据库基序。需要更多数据的示例（基因组、峰、reads 文件、背景集合）用带固定种子的随机数生成器生成，并写入临时目录。没有任何内容会写到该目录之外。

**示例文件：** `motifs.pfm`

```text
>MA0001.1 AR1
0.2 0.2 0.2 0.4
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
>MA0002.1 GATA
0.1 0.1 0.7 0.1
0.1 0.7 0.1 0.1
0.7 0.1 0.1 0.1
0.1 0.1 0.1 0.7
0.4 0.2 0.2 0.2
0.2 0.4 0.2 0.2
>MA0001.2 AR1b
0.2 0.2 0.2 0.4
0.1 0.1 0.7 0.1
0.7 0.1 0.1 0.1
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.6 0.2 0.1
```

**示例文件：** `seqs.fa`

```text
>s1
ACGTACGTAATTCCGGTGACACAGTCCATGACTTAGC
>s2
GGATCCATGACCTGACACTTAGGCATACCGATCAGTA
>s3
TTAGCGCATACGGATCCTGACACGTACTGAGCATCGA
>s4
CCGATTAGCGTACGTTAGCCATGCATGGATCAGTTCA
>s5
AGCATACTTGACACGGTACCGATAGCTTAGGCATCCA
>s6
TACCGTAGCTAGGATCCGATAGCTTACGGATCCGTAA
```

**示例文件：** `db.pfm`

```text
>DB1 long
0.25 0.25 0.25 0.25
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
0.25 0.25 0.25 0.25
0.25 0.25 0.25 0.25
```

## 1) `dnakit.sequence.gimmemotifs_cluster_motif_report()`

- **作用：** 按相似度对基序聚类，并写出带序列标志的 HTML 报告。
- **处理方法：** 创建输出目录（含 `images` 子目录），然后用相似度阈值对 PFM 文件调用 GimmeMotifs 的 `cluster_motifs_with_report()`。聚类后的 PFM 写入该目录下的 `clustered_motifs.pfm`。
- **API：** `gimmemotifs_cluster_motif_report(pfm_file: 'str', *, output_dir: 'str', threshold: 'float', title: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。
- **输入：** `pfm_file` 是已存在的 PFM 文件。`output_dir` 不能已存在。`threshold` 是 (0, 1] 内的数。`title` 是可选的报告标题。
- **输出：** `clusters`（列表，每项为含基序 ID 的 `{"representative", "members"}`）、`cluster_count`，以及 `outputs`（含 `clustered_pfm`，GimmeMotifs 写出报告时还有 `report`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_cluster_motif_report

out = Path(tempfile.mkdtemp()) / "report"
result = gimmemotifs_cluster_motif_report(
    "motifs.pfm", output_dir=out, threshold=0.9, verbose=False
)
print(result["cluster_count"], sorted(result["outputs"]))
for cluster in sorted(result["clusters"], key=lambda c: c["representative"]):
    print(cluster["representative"], sorted(cluster["members"]))
```

- **示例输出：**

```text
2 ['clustered_pfm', 'report']
Average_1 ['MA0001.1 AR1', 'MA0001.2 AR1b']
MA0002.1 GATA ['MA0002.1 GATA']
```

## 2) `dnakit.sequence.gimmemotifs_cluster_pfm()`

- **作用：** 按相似度对 PFM 基序聚类，写出聚类后的 PFM 和 HTML 报告。
- **处理方法：** 用固定参数列表且不经 shell 运行 `gimme cluster` 命令。`gimme` 创建新的输出目录，写出 `clustered_motifs.pfm`、`gimme.clustered.html` 以及每个聚类的序列标志；函数会检查前两个文件是否存在。
- **API：** `gimmemotifs_cluster_pfm(input_pfm: 'str', *, output_dir: 'str', threshold: 'float | None' = None, compare_reverse_complement: 'bool' = True, threads: 'int' = 1, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）及其 `gimme` 程序（依次从 `gimme_path`、PATH、当前 Python 解释器所在目录查找）。
- **输入：** `input_pfm` 是已存在的 PFM 文件。`output_dir` 不能已存在。`threshold` 为 None（用 GimmeMotifs 默认值）或 (0, 1] 内的数。`compare_reverse_complement=False` 会传入 `-s`（不比较反向互补）。`threads` 是正整数。
- **输出：** `outputs`（含 `clustered_pfm` 和 `report`）以及 `log_tail`（程序输出的最后 2000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_cluster_pfm

out = Path(tempfile.mkdtemp()) / "clusters"
result = gimmemotifs_cluster_pfm("motifs.pfm", output_dir=out, threshold=0.9, verbose=False)
print(sorted(result["outputs"]))
print((out / "clustered_motifs.pfm").read_text().count(">"))
```

- **示例输出：**

```text
['clustered_pfm', 'report']
2
```

## 3) `dnakit.sequence.gimmemotifs_combine_peak_summits()`

- **作用：** 把 MACS2 的 summit 或 narrowPeak 文件合并成一组互不重叠的 summit。
- **处理方法：** 用窗口宽度对这些文件调用 GimmeMotifs 的 `combine_peaks()`。每个窗口内只保留得分最高的 summit。GimmeMotifs 会为基因组 FASTA 建立索引，索引文件（`.fai`、`.sizes` 和 `.gaps.bed`）写在它旁边。行会被返回，也可以选择写成无表头的制表符分隔文件。
- **API：** `gimmemotifs_combine_peak_summits(peak_files: 'Sequence[str]', *, genome_fasta: 'str', window: 'int' = 200, scale_value: 'bool' = False, output_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。
- **输入：** `peak_files` 是非空的已存在 summit 或 narrowPeak 文件列表。`genome_fasta` 是已存在的 FASTA 文件。`window` 是正整数。`scale_value` 会按文件重新缩放得分。给出 `output_path` 时它不能已存在。
- **输出：** `columns`、`rows`（每个保留的 summit 一个字典）、`row_count`，给出 `output_path` 时还有含 `table` 的 `outputs`。`rows` 中从文件读出的文本列的值是字符串。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_combine_peak_summits

t = Path(tempfile.mkdtemp())
rng = random.Random(1)
genome = "".join(
    f">{c}\n{''.join(rng.choice('ACGT') for _ in range(3000))}\n" for c in ("chr1", "chr2")
)
(t / "genome.fa").write_text(genome)
(t / "a.bed").write_text("chr1\t500\t501\ta1\t12.5\nchr1\t1500\t1501\ta2\t8.0\nchr2\t700\t701\ta3\t5.0\n")
(t / "b.bed").write_text("chr1\t520\t521\tb1\t20.0\nchr2\t2000\t2001\tb2\t9.0\n")
result = gimmemotifs_combine_peak_summits(
    [str(t / "a.bed"), str(t / "b.bed")], genome_fasta=str(t / "genome.fa"), verbose=False
)
print(result["columns"], result["row_count"])
for row in result["rows"]:
    print(row["chrom"], row["start"], row["end"])
```

- **示例输出：**

```text
['chrom', 'start', 'end'] 4
chr1 520 521
chr1 1500 1501
chr2 700 701
chr2 2000 2001
```

## 4) `dnakit.sequence.gimmemotifs_coverage_table()`

- **作用：** 统计 BAM 或 BED 文件的 reads 在每个峰区域内的数量。
- **处理方法：** 先把峰复制到带 `+` 链列的临时 BED 文件，因为上游对其他链值会反转窗口并随后报错；只使用峰中心，所以结果不变。然后用窗口、对数变换、归一化、选取前若干区域、去重和去重复序列的设置调用 GimmeMotifs 的 `coverage_table()`。BAM 文件必须已有 `.bai` 索引，因为上游的自动建索引引用了未定义的名称。
- **API：** `gimmemotifs_coverage_table(peak_file: 'str', data_files: 'Sequence[str]', *, output_path: 'str | None' = None, window: 'int' = 200, log_transform: 'bool' = True, normalization: 'str | None' = None, top: 'int' = -1, top_method: 'str' = 'var', remove_duplicates: 'bool' = True, remove_repeats: 'bool' = True, threads: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它依赖 pysam 和 pybedtools。
- **输入：** `peak_file` 是至少含染色体、起点、终点的 BED 文件。`data_files` 是非空的 BAM（需索引）或 BED 文件列表。`normalization` 为 None、`"quantile"` 或 `"scale"`。`top_method` 为 `"var"`、`"std"`、`"mean"` 或 `"random"`；`top=-1` 保留所有区域。给出 `output_path` 时它不能已存在。
- **输出：** `columns`（每个数据文件一列）、`rows`（每个区域一个字典，含形如 `chr1:400-600` 的 `index` 以及每个文件一个值）、`row_count`，给出 `output_path` 时还有含 `table` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_coverage_table

t = Path(tempfile.mkdtemp())
(t / "peaks.bed").write_text("chr1\t400\t600\tp1\t1\nchr1\t1400\t1600\tp2\t1\nchr2\t600\t800\tp3\t1\n")
reads1 = "".join(f"chr1\t{450 + i * 10}\t{480 + i * 10}\tr{i}\t0\t+\n" for i in range(8))
reads2 = "".join(f"chr1\t{1450 + i * 10}\t{1480 + i * 10}\tr{i}\t0\t+\n" for i in range(5))
(t / "reads1.bed").write_text(reads1 + "chr2\t650\t680\tq1\t0\t+\n")
(t / "reads2.bed").write_text(reads2 + "chr2\t650\t680\tq1\t0\t+\nchr2\t660\t690\tq2\t0\t+\n")
result = gimmemotifs_coverage_table(
    str(t / "peaks.bed"),
    [str(t / "reads1.bed"), str(t / "reads2.bed")],
    log_transform=False,
    remove_duplicates=False,
    remove_repeats=False,
    verbose=False,
)
print(result["columns"], result["row_count"])
for row in sorted(result["rows"], key=lambda r: r["index"]):
    print(row["index"], row["reads1"], row["reads2"])
```

- **示例输出：**

```text
['reads1', 'reads2'] 3
chr1:1400-1600 0.0 5.0
chr1:400-600 8.0 0.0
chr2:600-800 1.0 2.0
```

## 5) `dnakit.sequence.gimmemotifs_create_background()`

- **作用：** 为基序分析创建一组背景序列（或区域）。
- **处理方法：** 先检查参数，因为上游遇到无效参数会直接退出进程，然后调用 GimmeMotifs 的 `create_background_file()`。`random` 基于 `input_file` 建立一阶马尔可夫模型；`genomic` 从基因组随机取长度为 `size`（或 `input_file` 序列长度中位数）的区域；`gc` 取与 `input_file` GC 含量匹配的基因组区域，需要 `bedtools` 程序来计算 GC 窗口。`promoter` 类型需要 genomepy 注释，因此被拒绝。
- **API：** `gimmemotifs_create_background(output_path: 'str', *, background_type: 'str', fmt: 'str' = 'fasta', size: 'int | None' = None, genome: 'str | None' = None, input_file: 'str | None' = None, number: 'int' = 10000, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。`gc` 类型还需要 PATH 中有 `bedtools`。
- **输入：** `output_path` 不能已存在。`background_type` 为 `random`、`genomic` 或 `gc`。`fmt` 为 `fasta` 或 `bed`（随机背景只能是 FASTA）。`random` 需要 `input_file`；`genomic` 和 `gc` 需要 `genome`（FASTA 路径）；`genomic` 还需要 `size` 或 `input_file`；`gc` 还需要 `input_file`。`number` 是正整数。
- **输出：** `outputs`（含 `background`）和 `sequence_count`，后者就是请求的 `number`（函数不会重新统计文件）。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_create_background

t = Path(tempfile.mkdtemp())
rng = random.Random(1)
genome = "".join(
    f">{c}\n{''.join(rng.choice('ACGT') for _ in range(3000))}\n" for c in ("chr1", "chr2")
)
(t / "genome.fa").write_text(genome)
result = gimmemotifs_create_background(
    t / "background.fa",
    background_type="genomic",
    genome=str(t / "genome.fa"),
    size=100,
    number=20,
    verbose=False,
)
records = (t / "background.fa").read_text().split(">")[1:]
print(result["sequence_count"], len(records), {len(r.split("\n", 1)[1].strip()) for r in records})
```

- **示例输出：**

```text
20 20 {100}
```

## 6) `dnakit.sequence.gimmemotifs_denovo_motifs()`

- **作用：** 用多个外部工具组成的集成流程，从 FASTA、BED 或 narrowPeak 文件预测从头（de novo）基序，再过滤、聚类并生成报告。
- **处理方法：** 检查输入文件存在，并检查至少有一个 GimmeMotifs 已知的从头预测程序已配置且可执行，否则在调用 GimmeMotifs 之前抛出 `BackendUnavailableError`。然后用输出目录、参数字典和三个开关调用 GimmeMotifs 的 `gimme_motifs()`。输出目录下写出的每个文件都按 `file_0`、`file_1` 等列出。
- **API：** `gimmemotifs_denovo_motifs(input_file: 'str', *, output_dir: 'str', params: 'Mapping[str, Any] | None' = None, filter_significant: 'bool' = True, cluster: 'bool' = True, create_report: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并至少有一个已配置的从头预测程序（例如 PATH 中的 MEME 或 HOMER，并登记在 `~/.config/gimmemotifs/gimmemotifs.cfg`）。默认参数要求 `gc` 背景，它需要基因组（`params["genome"]`）；没有基因组时 GimmeMotifs 会记录 “Need a genome to create background” 并让整个进程退出。完整的从头预测需要这里无法假设的工具和基因组，所以下面的示例演示输入检查：缺失的输入文件会被 `ConfigurationError` 拒绝。
- **输入：** `input_file` 是已存在的文件。`output_dir` 不能已存在。`params` 遵循 GimmeMotifs 的参数字典（`tools`、`analysis`、`background`、`genome`、`fraction` 等）。
- **输出：** `outputs`，按相对路径排序列出所有写出的文件。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_denovo_motifs

t = Path(tempfile.mkdtemp())
try:
    gimmemotifs_denovo_motifs(str(t / "missing.fa"), output_dir=t / "denovo", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(": ")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_file does not exist
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_denovo_motifs

result = gimmemotifs_denovo_motifs(
    "peaks.fa",
    output_dir="denovo_out",
    params={"tools": "MEME,Homer", "background": "random", "analysis": "small"},
)
print(result["outputs"].keys())
```

## 7) `dnakit.sequence.gimmemotifs_denovo_predict()`

- **作用：** 在 FASTA 文件上对照背景 FASTA 运行已配置的从头预测工具，并写出预测到的基序。
- **处理方法：** 检查两个输入文件和新的输出路径，检查至少有一个从头预测程序已配置且可执行，然后调用 GimmeMotifs 的 `predict_motifs()`，它会运行 `params["tools"]` 中指定的工具。若没有写出基序文件，函数抛出 `BackendExecutionError`。
- **API：** `gimmemotifs_denovo_predict(input_fasta: 'str', background_fasta: 'str', output_path: 'str', *, params: 'Mapping[str, Any] | None' = None, stats_fg: 'Mapping[str, Any] | None' = None, stats_bg: 'Mapping[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并至少有一个已配置的从头预测程序。这里无法假设这些工具，所以下面的示例演示输入检查：缺失的背景文件会被 `ConfigurationError` 拒绝。
- **输入：** `input_fasta` 和 `background_fasta` 是已存在的 FASTA 文件。`output_path` 不能已存在。`params` 是 GimmeMotifs 的参数字典；`stats_fg` 和 `stats_bg` 是可选的统计输入，原样传入。
- **输出：** `outputs`（含 `motifs`，即 GimmeMotifs 写出的 PFM 文件路径）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_denovo_predict

t = Path(tempfile.mkdtemp())
try:
    gimmemotifs_denovo_predict("seqs.fa", str(t / "missing.fa"), t / "motifs.pfm", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(": ")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] background_fasta does not exist
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_denovo_predict

result = gimmemotifs_denovo_predict(
    "peaks.fa",
    "background.fa",
    "predicted.pfm",
    params={"tools": "MEME,Homer", "analysis": "small", "genome": "hg38"},
)
print(result["outputs"]["motifs"])
```

## 8) `dnakit.sequence.gimmemotifs_gc_bin_bed()`

- **作用：** 在基因组的各个 GC 含量区间内均匀抽取 BED 区域。
- **处理方法：** 调用 GimmeMotifs 的 `gc_bin_bedfile()`。它用 pybedtools 计算基因组窗口的 GC 含量，并把该索引缓存到 `<缓存目录>/gimmemotifs/<基因组文件名>.gcfreq.<min_bin_size>.feather`。缓存以文件名而不是路径为键，所以同名的不同基因组会复用旧索引。缓存目录由 `XDG_CACHE_HOME` 决定；示例把它设为临时目录。
- **API：** `gimmemotifs_gc_bin_bed(output_path: 'str', *, genome: 'str', number: 'int', length: 'int' = 200, bins: 'Sequence[float] | None' = None, random_state: 'int | None' = None, min_bin_size: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并且 PATH 中有 `bedtools` 程序（通过 pybedtools 调用）。
- **输入：** `output_path` 不能已存在。`genome` 是已存在的 FASTA 文件。`number` 是正整数，且不小于区间数；每个区间有放回地抽取 `number // 区间数` 个区域。`length` 是区域长度。`bins` 是 `(low, high)` GC 范围的列表；默认是 0-0.2、以 0.05 为步长的 0.2-0.8，以及 0.8-1。`random_state` 使抽样可重复。
- **输出：** `outputs`（含 `bed`）和 `region_count`，后者是请求的 `number`。BED 文件含染色体、起点、终点和区间标签。
- **示例代码：**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_gc_bin_bed

rng = random.Random(1)
with (t / "genome.fa").open("w") as handle:
    for chrom in ("chr1", "chr2"):
        seq = ""
        for block in range(24):
            gc = (0.3, 0.4, 0.5, 0.6)[block % 4]
            seq += "".join(
                rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(500)
            )
        handle.write(f">{chrom}\n{seq}\n")
result = gimmemotifs_gc_bin_bed(
    t / "bins.bed",
    genome=str(t / "genome.fa"),
    number=12,
    length=100,
    bins=[(0.2, 0.4), (0.4, 0.6)],
    random_state=1,
    verbose=False,
)
rows = [line.split("\t") for line in (t / "bins.bed").read_text().splitlines()]
print(result["region_count"], len(rows), sorted({row[3] for row in rows}))
```

- **示例输出：**

```text
12 12 ['0.20-0.40', '0.40-0.60']
```

## 9) `dnakit.sequence.gimmemotifs_gc_bin_index()`

- **作用：** 写出固定大小基因组窗口的 GC 含量索引，用于 GC 匹配抽样。
- **处理方法：** 调用 GimmeMotifs 的 `create_gc_bin_index()`。它用 pybedtools 把基因组切成 `min_bin_size` 大小的窗口，计算各窗口的 GC 含量和 N 数量，再加上 2 个和 5 个窗口的平均值，并写成 Feather 文件。函数会检查文件已创建。
- **API：** `gimmemotifs_gc_bin_index(genome: 'str', output_path: 'str', *, min_bin_size: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并且 PATH 中有 `bedtools` 程序（通过 pybedtools 调用）。
- **输入：** `genome` 是已存在的 FASTA 文件。`output_path` 不能已存在。`min_bin_size` 是正整数；取值很小会使索引非常大。
- **输出：** `outputs`（含 `index`，即 Feather 文件）。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

import pandas

from dnakit.sequence import gimmemotifs_gc_bin_index

t = Path(tempfile.mkdtemp())
rng = random.Random(1)
(t / "genome.fa").write_text(f">chr1\n{''.join(rng.choice('ACGT') for _ in range(2000))}\n")
result = gimmemotifs_gc_bin_index(str(t / "genome.fa"), t / "gc.feather", verbose=False)
table = pandas.read_feather(t / "gc.feather")
print(sorted(result["outputs"]), len(table), list(table.columns))
```

- **示例输出：**

```text
['index'] 20 ['chrom', 'start', 'end', 'w100', 'n100', 'w200', 'n200', 'w500', 'n500']
```

## 10) `dnakit.sequence.gimmemotifs_histogram_plot()`

- **作用：** 把一组数字画成 SVG 直方图。
- **处理方法：** 把值转成浮点数，输出路径缺少 `.svg` 时自动补上，然后用分箱数、范围、标题和横轴标签调用 GimmeMotifs 的 `plot_histogram()`。
- **API：** `gimmemotifs_histogram_plot(values: 'Sequence[float]', output_path: 'str', *, x_range: 'tuple[float, float] | None' = None, bins: 'int' = 10, title: 'str | None' = None, x_label: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib。
- **输入：** `values` 是非空的数字列表。补上 `.svg` 之后 `output_path` 不能已存在。`x_range` 为 None 或满足 `low < high` 的 `(low, high)`。`bins` 是正整数。
- **输出：** `outputs`（含 `image`，即实际写出的 SVG 路径）、`count`，以及这些值的 `minimum` 和 `maximum`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_histogram_plot

t = Path(tempfile.mkdtemp())
result = gimmemotifs_histogram_plot(
    [0.1, 0.2, 0.2, 0.3, 0.5, 0.9], t / "hist", bins=5, title="scores", verbose=False
)
image = Path(result["outputs"]["image"])
print(result["count"], result["minimum"], result["maximum"], image.name, image.read_text()[:4])
```

- **示例输出：**

```text
6 0.1 0.9 hist.svg <?xm
```

## 11) `dnakit.sequence.gimmemotifs_jaspar_motif_info()`

- **作用：** 从 JASPAR 网络服务获取某个 JASPAR 矩阵的元数据。
- **处理方法：** 检查 ID 是非空字符串，然后调用 GimmeMotifs 的 `get_jaspar_motif_info()`，它请求 `http://jaspar.genereg.net/api/v1/matrix/<id>?format=json` 并返回 JSON 对象。这是本组中唯一使用网络的函数。
- **API：** `gimmemotifs_jaspar_motif_info(motif_id: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并且能访问互联网上的 `jaspar.genereg.net`；不需要密钥。下面的示例不访问网络，所以演示输入检查：空 ID 会被 `ConfigurationError` 拒绝。
- **输入：** `motif_id` 是 JASPAR 矩阵 ID，例如 `MA0106.3`，必须是非空字符串。
- **输出：** `motif_id` 和 `info`（JASPAR 返回的 JSON 对象）。
- **示例代码：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_jaspar_motif_info

try:
    gimmemotifs_jaspar_motif_info("  ", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error))
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] motif_id must be a non-empty string.
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_jaspar_motif_info

result = gimmemotifs_jaspar_motif_info("MA0106.3")
print(result["info"]["name"])
```

## 12) `dnakit.sequence.gimmemotifs_logo_images()`

- **作用：** 为 PFM 文件中的基序画序列标志，保存为 PNG 文件。
- **处理方法：** 在新的输出目录中运行 `gimme logo` 命令：用 `-i` 传入所选 ID，用 `-k` 传入标志类型，关闭标题时加 `--notitle`。返回目录下找到的每个 PNG；一个都没有时函数抛出 `BackendExecutionError`。
- **API：** `gimmemotifs_logo_images(pfm_file: 'str', *, output_dir: 'str', motif_ids: 'Sequence[str] | None' = None, kind: 'str | None' = None, show_title: 'bool' = True, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）及其 `gimme` 程序（依次从 `gimme_path`、PATH、当前 Python 解释器所在目录查找）。
- **输入：** `pfm_file` 是已存在的 PFM 文件。`output_dir` 不能已存在。`motif_ids` 选择基序（默认全部）。`kind` 为 None、`information`、`frequency`、`energy` 或 `ensembl`。
- **输出：** `outputs`（`image_0`、`image_1` 等，为排序后的 PNG 路径）以及 `log_tail`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_logo_images

out = Path(tempfile.mkdtemp()) / "logos"
result = gimmemotifs_logo_images(
    "motifs.pfm", output_dir=out, motif_ids=["MA0001.1 AR1", "MA0002.1 GATA"], verbose=False
)
print(sorted(Path(path).name for path in result["outputs"].values()))
```

- **示例输出：**

```text
['MA0001.1 AR1.png', 'MA0002.1 GATA.png']
```

## 13) `dnakit.sequence.gimmemotifs_match_motifs()`

- **作用：** 为每个查询基序找出数据库中最匹配的基序。
- **处理方法：** 运行 `gimme match -d <数据库> -n <nreport> <查询>`，把它的表格流式写入新的输出文件，再把表格解析成行。不提供 `gimme match -o` 的图形报告，因为它在当前 GimmeMotifs 版本里会因 numpy 广播错误而失败。
- **API：** `gimmemotifs_match_motifs(pfm_file: 'str', *, database_pfm: 'str', output_path: 'str', nreport: 'int' = 1, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）及其 `gimme` 程序（依次从 `gimme_path`、PATH、当前 Python 解释器所在目录查找）。
- **输入：** `pfm_file` 是查询 PFM 文件，`database_pfm` 是数据库 PFM 文件，二者都必须存在。`output_path` 不能已存在。`nreport` 是每个查询返回的匹配数（正整数）。
- **输出：** `rows`（每个匹配一个字典，含表格列 `Motif`、`Match`、`Score`、`P-value`，值均为字符串）、`row_count`、含 `table` 的 `outputs`，以及 `log_tail`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_match_motifs

out = Path(tempfile.mkdtemp()) / "match.tsv"
result = gimmemotifs_match_motifs("motifs.pfm", database_pfm="db.pfm", output_path=out, verbose=False)
print(result["row_count"])
for row in result["rows"]:
    print(row["Motif"], "->", row["Match"], row["Score"])
```

- **示例输出：**

```text
3
MA0001.1 AR1 -> DB1 long 0.99
MA0002.1 GATA -> DB1 long 0.23
MA0001.2 AR1b -> DB1 long 0.96
```

## 14) `dnakit.sequence.gimmemotifs_match_plot()`

- **作用：** 把每个查询基序与其最佳数据库匹配及 p 值并排画出。
- **处理方法：** 按 ID 从各自的 PFM 文件读出每个查询基序和数据库基序，然后用 `(查询, 匹配, p 值)` 三元组列表调用 GimmeMotifs 的 `match_plot()`，把图写到新的输出路径。
- **API：** `gimmemotifs_match_plot(matches: 'Sequence[Mapping[str, Any]]', output_path: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib。
- **输入：** `matches` 是非空的字典列表，每项含 `query_pfm`、`query_id`、`database_pfm`、`database_id`（已存在的文件及其中的基序 ID）和数值型 `p_value`。`output_path` 不能已存在。
- **输出：** `outputs`（含 `image`）和 `pairs`（绘制的配对数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_match_plot

out = Path(tempfile.mkdtemp()) / "match.png"
result = gimmemotifs_match_plot(
    [
        {
            "query_pfm": "motifs.pfm",
            "query_id": "MA0001.1 AR1",
            "database_pfm": "db.pfm",
            "database_id": "DB1 long",
            "p_value": 0.01,
        }
    ],
    out,
    verbose=False,
)
print(result["pairs"], out.read_bytes()[:4])
```

- **示例输出：**

```text
1 b'\x89PNG'
```

## 15) `dnakit.sequence.gimmemotifs_matched_gc_bed()`

- **作用：** 抽取 GC 含量与另一个 BED 文件区域相仿的 BED 区域。
- **处理方法：** 调用 GimmeMotifs 的 `matched_gc_bedfile()`。它读取 `match_bed`（BED 或 FASTA 文件）中区域的 GC 含量，从基因组的 GC 索引中抽取窗口，并写出抽到的区域。索引按基因组文件名缓存在缓存目录中，见第 8 节。
- **API：** `gimmemotifs_matched_gc_bed(output_path: 'str', *, match_bed: 'str', genome: 'str', number: 'int', size: 'int | None' = None, min_bin_size: 'int' = 100, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并且 PATH 中有 `bedtools` 程序（通过 pybedtools 调用）。
- **输入：** `output_path` 不能已存在。`match_bed` 和 `genome` 是已存在的文件，基因组是 FASTA 文件。`number` 是正整数。`size` 是区域长度（默认用输入区域的长度）。`min_bin_size` 是正整数。
- **输出：** `outputs`（含 `bed`）。函数不会统计区域数。
- **示例代码：**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_matched_gc_bed

rng = random.Random(1)
with (t / "genome.fa").open("w") as handle:
    for chrom in ("chr1", "chr2"):
        seq = ""
        for block in range(24):
            gc = (0.3, 0.4, 0.5, 0.6)[block % 4]
            seq += "".join(
                rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(500)
            )
        handle.write(f">{chrom}\n{seq}\n")
(t / "match.bed").write_text("chr1\t0\t100\nchr1\t500\t600\nchr2\t1000\t1100\n")
result = gimmemotifs_matched_gc_bed(
    t / "matched.bed",
    match_bed=str(t / "match.bed"),
    genome=str(t / "genome.fa"),
    number=10,
    verbose=False,
)
rows = [line.split("\t") for line in (t / "matched.bed").read_text().splitlines()]
print(sorted(result["outputs"]), len(rows), {int(r[2]) - int(r[1]) for r in rows})
```

- **示例输出：**

```text
['bed'] 10 {100}
```

## 16) `dnakit.sequence.gimmemotifs_motif_image_series()`

- **作用：** 为列出的每个基序 ID 写出一张序列标志 PNG。
- **处理方法：** 读取 PFM 文件中的基序，检查每个 ID 都存在，然后调用 GimmeMotifs 的 `motif_to_img_series()`，在 `<output_dir>/<subdir>` 下为每个基序画一张 PNG。ID 中除字母、数字和 `-` 以外的字符在文件名里被替换为 `_`。函数会检查每张图片都已写出。
- **API：** `gimmemotifs_motif_image_series(motif_ids: 'Sequence[str]', *, pfm_file: 'str', output_dir: 'str', subdir: 'str' = 'logos', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），以及 pandas。
- **输入：** `motif_ids` 是非空的 ID 列表，这些 ID 必须在 `pfm_file` 中。`output_dir` 不能已存在。`subdir` 是存放图片的子目录。
- **输出：** `images`（基序 ID 对应 PNG 路径）和 `image_count`。`outputs` 为空。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_motif_image_series

out = Path(tempfile.mkdtemp()) / "images"
result = gimmemotifs_motif_image_series(
    ["MA0001.1 AR1", "MA0002.1 GATA"], pfm_file="motifs.pfm", output_dir=out, verbose=False
)
print(result["image_count"])
for motif_id, path in result["images"].items():
    print(motif_id, "->", Path(path).relative_to(out))
```

- **示例输出：**

```text
2
MA0001.1 AR1 -> logos/MA0001_1_AR1.png
MA0002.1 GATA -> logos/MA0002_1_GATA.png
```

## 17) `dnakit.sequence.gimmemotifs_motif_localization()`

- **作用：** 画出基序在等长序列上的出现位置，并检验这些位置是否均匀分布。
- **处理方法：** 按 ID 读取基序并调用 GimmeMotifs 的 `motif_localization()`：它扫描 FASTA 文件，把匹配位置画成 SVG 直方图，并返回相对均匀分布的 Kolmogorov-Smirnov p 值。上游只在基序有匹配时才写 SVG；否则 `outputs` 中没有 `image`。
- **API：** `gimmemotifs_motif_localization(fasta_file: 'str', *, pfm_file: 'str', motif_id: 'str', output_path: 'str', size: 'int', cutoff: 'float' = 0.9, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib。
- **输入：** `fasta_file` 和 `pfm_file` 是已存在的文件；`motif_id` 是 PFM 文件中的 ID。`output_path` 不能已存在。`size` 是直方图使用的序列长度，必须大于基序长度。`cutoff` 是得分阈值，取 (0, 1] 内的数。
- **输出：** `motif_id`、`p_value`，以及写出 SVG 时 `outputs` 中的 `image`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_motif_localization

out = Path(tempfile.mkdtemp()) / "location.svg"
result = gimmemotifs_motif_localization(
    "seqs.fa",
    pfm_file="motifs.pfm",
    motif_id="MA0001.1 AR1",
    output_path=out,
    size=37,
    cutoff=0.9,
    verbose=False,
)
print(result["motif_id"], round(result["p_value"], 3), sorted(result["outputs"]), out.exists())
```

- **示例输出：**

```text
MA0001.1 AR1 0.244 ['image'] True
```

## 18) `dnakit.sequence.gimmemotifs_motif_location_plots()`

- **作用：** 画出每个基序相对于序列中心的出现位置。
- **处理方法：** 在新的输出目录中运行 `gimme location`，传入得分阈值以及可选的大小和基序 ID。返回目录下找到的每个文件。
- **API：** `gimmemotifs_motif_location_plots(fasta_file: 'str', *, pfm_file: 'str', output_dir: 'str', size: 'int | None' = None, motif_ids: 'Sequence[str] | None' = None, cutoff: 'float' = 0.95, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）及其 `gimme` 程序（依次从 `gimme_path`、PATH、当前 Python 解释器所在目录查找）。
- **输入：** `fasta_file` 和 `pfm_file` 是已存在的文件。`output_dir` 不能已存在。`size` 是正整数（绘图所用的序列长度）。`motif_ids` 选择基序。`cutoff` 是 (0, 1] 内的数。
- **输出：** `outputs`（`file_0`、`file_1` 等，为排序后的路径；每个基序一张直方图 SVG）以及 `log_tail`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_motif_location_plots

out = Path(tempfile.mkdtemp()) / "location"
result = gimmemotifs_motif_location_plots(
    "seqs.fa", pfm_file="motifs.pfm", output_dir=out, size=36, verbose=False
)
print(sorted(Path(path).name for path in result["outputs"].values()))
```

- **示例输出：**

```text
['MA0001.1 AR1_histogram.svg', 'MA0001.2 AR1b_histogram.svg', 'MA0002.1 GATA_histogram.svg']
```

## 19) `dnakit.sequence.gimmemotifs_parallel_denovo_predict()`

- **作用：** 用多个外部从头预测工具并行预测基序。
- **处理方法：** 检查 FASTA 文件和新的输出路径，检查至少有一个从头预测程序已配置且可执行，然后用分析规模、物种、工具列表、线程数和时间限制调用 GimmeMotifs 的 `pp_predict_motifs()`。若没有写出基序文件，函数抛出 `BackendExecutionError`。
- **API：** `gimmemotifs_parallel_denovo_predict(fasta_file: 'str', output_path: 'str', *, analysis: 'str' = 'small', organism: 'str' = 'hg19', single: 'bool' = False, background: 'str' = '', tools: 'str | None' = None, threads: 'int' = 8, max_time: 'int' = -1, stats_fg: 'Mapping[str, Any] | None' = None, stats_bg: 'Mapping[str, Any] | None' = None, gc: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并至少有一个已配置的从头预测程序。下面的示例不运行任何工具：缺失的输入文件在任何程序启动之前就会被 `ConfigurationError` 拒绝。
- **输入：** `fasta_file` 是已存在的 FASTA 文件。`output_path` 不能已存在。`analysis` 是 GimmeMotifs 的分析规模（`xs`、`small`、`medium`、`large`、`xl`）。`tools` 是以逗号分隔的工具名（为 None 时用 GimmeMotifs 默认集合）。`threads` 是正整数。`max_time` 是以秒计的时间限制（-1 表示不限）。
- **输出：** `outputs`（含 `motifs`，即 GimmeMotifs 写出的 PFM 文件）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_parallel_denovo_predict

t = Path(tempfile.mkdtemp())
try:
    gimmemotifs_parallel_denovo_predict(str(t / "missing.fa"), t / "motifs.pfm", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(": ")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] fasta_file does not exist
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_parallel_denovo_predict

result = gimmemotifs_parallel_denovo_predict(
    "peaks.fa", "predicted.pfm", analysis="small", tools="MEME,Homer", threads=4
)
print(result["outputs"]["motifs"])
```

## 20) `dnakit.sequence.gimmemotifs_pfm_scan()`

- **作用：** 用固定得分阈值，以 PFM 基序扫描 FASTA 或 BED 序列。
- **处理方法：** 运行 `gimme scan -p <pfm> -c <阈值> -n <nreport> -N <线程数>`，并按格式加 `-t`（表格）、`-b`（BED）或 `-T`（得分表），关闭反向互补时再加 `-r`，把结果流式写入新输出目录中的 `scan.tsv` 或 `scan.bed`。阈值是必需的，因为 `gimme scan` 的 FPR 模式需要背景集或基因组，而本封装不接受它们。表格输出会被解析成行。
- **API：** `gimmemotifs_pfm_scan(input_file: 'str', *, pfm_file: 'str', output_dir: 'str', cutoff: 'float', nreport: 'int' = 1, reverse_complement: 'bool' = True, output_format: 'str' = 'table', threads: 'int' = 1, gimme_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）及其 `gimme` 程序（依次从 `gimme_path`、PATH、当前 Python 解释器所在目录查找）。
- **输入：** `input_file` 是 FASTA 或 BED 文件，`pfm_file` 是 PFM 文件，二者都必须存在。`output_dir` 不能已存在。`cutoff` 是 (0, 1] 内的数（基序得分阈值）。`output_format` 为 `table`、`bed` 或 `score_table`。`nreport` 和 `threads` 是正整数。
- **输出：** `outputs`（含 `scan`）、`rows`（每条序列一个字典；表格格式下是每个基序的匹配数，字符串形式）、`row_count`（BED 输出时为 0）以及 `log_tail`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_pfm_scan

out = Path(tempfile.mkdtemp()) / "scan"
result = gimmemotifs_pfm_scan(
    "seqs.fa", pfm_file="motifs.pfm", output_dir=out, cutoff=0.9, verbose=False
)
print(result["row_count"], Path(result["outputs"]["scan"]).name)
for row in result["rows"]:
    print(row["sequence"], row["MA0001.1 AR1"], row["MA0002.1 GATA"], row["MA0001.2 AR1b"])
```

- **示例输出：**

```text
6 scan.tsv
s1 1 0 1
s2 1 1 1
s3 1 1 1
s4 0 0 0
s5 1 1 1
s6 0 0 0
```

## 21) `dnakit.sequence.gimmemotifs_roc_curve_plot()`

- **作用：** 画出一条或多条 ROC 曲线，横轴为 1 - 特异性。
- **处理方法：** 把输入转成浮点数组，路径没有图片扩展名时补上 `.png`，并以曲线 ID 作为图例调用 GimmeMotifs 的 `roc_plot()`。接受的扩展名是 GimmeMotifs 支持的那些（`.png`、`.pdf`、`.svg`、`.ps`）。
- **API：** `gimmemotifs_roc_curve_plot(output_path: 'str', *, plot_x: 'Sequence[Any]', plot_y: 'Sequence[Any]', curve_ids: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib。
- **输入：** `plot_x` 和 `plot_y` 非空且长度相同：可以是一条曲线（数字列表），也可以是多条曲线（列表的列表）。`curve_ids` 为曲线命名。补上扩展名之后 `output_path` 不能已存在。
- **输出：** `outputs`（含 `image`，即实际写出的路径）和 `curves`（曲线条数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_roc_curve_plot

out = Path(tempfile.mkdtemp()) / "roc"
result = gimmemotifs_roc_curve_plot(
    out,
    plot_x=[[0, 0.2, 1], [0, 0.4, 1]],
    plot_y=[[0, 0.8, 1], [0, 0.6, 1]],
    curve_ids=["motif A", "motif B"],
    verbose=False,
)
image = Path(result["outputs"]["image"])
print(result["curves"], image.name, image.read_bytes()[1:4])
```

- **示例输出：**

```text
2 roc.png b'PNG'
```

## 22) `dnakit.sequence.gimmemotifs_gimme_motifs()`

- **作用：** 对样本 FASTA 文件运行 `gimme motifs`：扫描已知基序（可选再做从头预测），对照背景计算富集统计，并写出报告。
- **处理方法：** 检查参数和从头预测工具（每个工具都必须已配置并已安装，否则在任何工作开始之前抛出 `BackendUnavailableError`），然后用显式构造的参数集调用 GimmeMotifs 的 `motifs` 命令函数。`pfm_file` 中的已知基序在样本和背景上被扫描；ROC、精确率-召回率和富集统计写入 `gimme.roc.report.txt`；选出一个非冗余子集，其扫描结果以 BED 文件写入 `motif_scan_results`；`report` 为真时写出 HTML 报告 `gimme.motifs.html` 和 `gimme.motifs.redundant.html`。基因组不是输入，所以 z 分数和 GC 归一化是关闭的。`background_fasta` 为 None 时，GimmeMotifs 从样本生成随机背景。在验证过的环境（xxhash 4.0.1）里，只要选出了至少一个非冗余基序，运行就会因 `TypeError: Strings must be encoded before hashing` 失败，因为后续的 BED 扫描用 `str` 对基序做哈希；示例使用的数据不会选出任何基序。
- **API：** `gimmemotifs_gimme_motifs(sample_fasta: 'str', pfm_file: 'str', output_dir: 'str', *, background_fasta: 'str | None' = None, known: 'bool' = True, denovo_tools: 'Sequence[str] | None' = None, analysis: 'str' = 'xl', fraction: 'float' = 0.2, single: 'bool' = False, report: 'bool' = True, ncpus: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。`denovo_tools` 中指定的从头预测工具必须已配置并已安装；默认（只用已知基序）不需要任何工具。
- **输入：** `sample_fasta`、`pfm_file` 和可选的 `background_fasta` 是已存在的文件。`output_dir` 不能已存在。`analysis` 为 `small`、`medium`、`large` 或 `xl`。`fraction` 是 (0, 1] 内的数。`known=False` 时必须给出 `denovo_tools`。`denovo_tools` 是工具名列表，例如 `["MEME", "Homer", "BioProspector"]`。`ncpus` 至少为 1。
- **输出：** `known`、`denovo_tools`、`motifs_scored`（ROC 报告中的基序数），以及 `outputs`（含 `roc_report`、`motif_scan_results`，存在时还有 `gimme.motifs.html`、`gimme.motifs.redundant.html` 和 `combined.pfm`）。
- **示例代码：**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_gimme_motifs

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


sample = [seq(200) for _ in range(40)]
sample = [s[:50] + "TGACGTCA" + s[58:] if i % 2 == 0 else s for i, s in enumerate(sample)]
(t / "sample.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(sample)))
(t / "background.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
result = gimmemotifs_gimme_motifs(
    str(t / "sample.fa"),
    "motifs.pfm",
    t / "run",
    background_fasta=str(t / "background.fa"),
    verbose=False,
)
print(result["known"], result["denovo_tools"], result["motifs_scored"])
print(sorted(result["outputs"]))
print((t / "run" / "gimme.roc.report.txt").read_text().splitlines()[0].split("\t")[:5])
```

- **示例输出：**

```text
True [] 3
['gimme.motifs.html', 'gimme.motifs.redundant.html', 'motif_scan_results', 'roc_report']
['Motif', '# matches', '% matches input', '# matches background', '% matches background']
```

## 23) `dnakit.sequence.gimmemotifs_maelstrom_html_report()`

- **作用：** 为已完成的 maelstrom 运行目录写出 HTML 报告。
- **处理方法：** 从运行目录读取 `final.out.txt`（基序对方法的表），统计在至少一个方法列中绝对值不低于阈值的基序数（列名含 `corr` 或 `% with motif` 的列不算方法列），然后调用 GimmeMotifs 的 `maelstrom_html_report()`。报告 `gimme.maelstrom.report.html` 写在运行目录里，该目录中已有的报告会被替换，与上游一致。
- **API：** `gimmemotifs_maelstrom_html_report(run_folder: 'str', pfm_file: 'str', *, threshold: 'float' = 3, plot_all_motifs: 'bool' = False, plot_no_motifs: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），以及 pandas。
- **输入：** `run_folder` 是含 `final.out.txt` 的已存在目录；`pfm_file` 是该次运行所用的基序数据库。`threshold` 不小于 0。`plot_all_motifs` 为每个基序画序列标志；`plot_no_motifs` 去掉标志列。
- **输出：** `threshold`、`motifs_in_report`、`methods`（方法列数），以及 `outputs`（含 `report`，目录存在时还有 `logos`）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_maelstrom_html_report

run = Path(tempfile.mkdtemp())
(run / "final.out.txt").write_text(
    "\tHypergeom\tMean\tHypergeom % with motif\tcorr Hypergeom\n"
    "MA0001.1 AR1\t5.2\t3.1\t80.0\t0.9\n"
    "MA0002.1 GATA\t-1.2\t0.2\t20.0\t0.1\n"
)
result = gimmemotifs_maelstrom_html_report(str(run), "motifs.pfm", threshold=3, verbose=False)
print(result["threshold"], result["motifs_in_report"], result["methods"], sorted(result["outputs"]))
html = (run / "gimme.maelstrom.report.html").read_text()
print("MA0001.1 AR1" in html, "MA0002.1 GATA" in html)
```

- **示例输出：**

```text
3.0 1 2 ['logos', 'report']
True False
```

## 24) `dnakit.sequence.gimmemotifs_predict_motifs()`

- **作用：** 对 FASTA 文件运行一个从头预测工具，并把预测到的基序写成 PFM 文件。
- **处理方法：** 在 GimmeMotifs 中查找该工具（它必须在 `gimmemotifs.cfg` 中配置且可执行），用 FASTA 文件和 `gimme prediction` 从参数文件读取的那种参数字典调用其 `run()` 方法，并用 `to_pfm()` 把每个预测到的基序依次写入新的输出文件。
- **API：** `gimmemotifs_predict_motifs(tool: 'str', fasta: 'str', output_pfm: 'str', *, params: 'Mapping[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并且所指定的工具已配置并已安装（例如 MEME 或 Homer；GimmeMotifs 也自带一些工具）。已知但未配置或不可执行的工具会抛出 `BackendUnavailableError`。下面的示例使用未知的工具名，它在任何程序启动之前就被拒绝。
- **输入：** `tool` 是 GimmeMotifs 的工具名（不区分大小写，须在其工具列表中）。`fasta` 是已存在的 FASTA 文件。`output_pfm` 不能已存在。`params` 是工具参数字典（例如 `background`、`width`、`number`、`analysis`、`use_strand`），也可以是 None。
- **输出：** `tool`、`motifs`（预测到的基序数）、`stdout_tail` 和 `stderr_tail`（工具输出的有限长度尾部），以及含 `pfm` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import gimmemotifs_predict_motifs

out = Path(tempfile.mkdtemp()) / "predicted.pfm"
try:
    gimmemotifs_predict_motifs("nosuchtool", "seqs.fa", out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(";")[0])
print(out.exists())
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] Unknown gimmemotifs tool 'nosuchtool'
False
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_predict_motifs

result = gimmemotifs_predict_motifs(
    "BioProspector",
    "peaks.fa",
    "predicted.pfm",
    params={"background": "background.fa", "width": 8, "number": 3, "analysis": "small", "use_strand": False},
)
print(result["motifs"])
```

## 25) `dnakit.sequence.gimmemotifs_annot2primpep()`

- **作用：** 写出一个基因组组装的主要肽段 FASTA：每个基因取最长的蛋白。
- **处理方法：** 检查 PATH 中有 `gffread`，把基因组 FASTA 和 GTF 注释复制到 `<output_dir>/<name>/`，然后调用 GimmeMotifs 的 `annot2primpep()`。它运行 gffread 翻译注释，并保留每个基因最长的蛋白。记录名为 `gene_name|gene_id`（基因名缺失时用 `.`）。结果是 `<output_dir>/prim_genes/<name>.pep.fa`。
- **API：** `gimmemotifs_annot2primpep(genome_fasta: 'str', annotation_gtf: 'str', output_dir: 'str', *, name: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并且 PATH 中有外部程序 `gffread`（例如 `conda install -c bioconda gffread`）。下面的示例在缺少 `gffread` 的情况下运行并展示报错；验证本页所用的环境没有安装该程序。
- **输入：** `genome_fasta` 和 `annotation_gtf` 是已存在的文件。`output_dir` 不能已存在。`name` 是不含 `/` 和 `'` 的普通标识符（默认用基因组文件名主干）。
- **输出：** `assembly`（名称）、`protein_count`（肽段文件的记录数）以及含 `peptides` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence import gimmemotifs_annot2primpep

out = Path(tempfile.mkdtemp()) / "assembly"
try:
    gimmemotifs_annot2primpep("seqs.fa", "motifs.pfm", out, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(" Hint")[0])
print(out.exists())
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] The gffread program was not found on PATH.
False
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_annot2primpep

result = gimmemotifs_annot2primpep("genome.fa", "genes.gtf", "primary_peptides", name="mygenome")
print(result["assembly"], result["protein_count"], result["outputs"]["peptides"])
```

## 26) `dnakit.sequence.gimmemotifs_best_motif_in_cluster()`

- **作用：** 返回每个从头预测基序聚类中的最佳基序。
- **处理方法：** 读取单个基序和聚类基序，解析每个聚类的 ID，对每个背景 FASTA 为所有相关基序计算统计（`calc_stats()`，关闭 z 分数和 GC 选项），再用这些统计调用 GimmeMotifs 的 `best_motif_in_cluster()`，它按给定指标对每个聚类内的基序排序。排序指标是秩统计量，因此不受逐基序 z 分数影响。
- **API：** `gimmemotifs_best_motif_in_cluster(single_motif_file: 'str', cluster_motif_file: 'str', clusters: 'Sequence[Mapping[str, Any]]', foreground_fasta: 'str', background_fastas: 'Mapping[str, str]', genome_fasta: 'str', *, metrics: 'Sequence[str]' = ('roc_auc', 'recall_at_fdr'), verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。
- **输入：** `single_motif_file` 和 `cluster_motif_file` 是 PFM 文件。`clusters` 是非空的字典列表，每项含 `representative`（聚类文件中的 ID）和 `members`（单个基序文件中的 ID）。`foreground_fasta`、`genome_fasta` 以及 `background_fastas`（背景名到 FASTA 路径的字典）中的每个路径都必须存在。`metrics` 是 GimmeMotifs 的统计名。
- **输出：** `best_motifs`，即每个聚类中最佳基序的 ID，按从好到差排列。
- **示例代码：**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_best_motif_in_cluster

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


foreground = [seq(200) for _ in range(40)]
foreground = [s[:50] + "TGACAC" + s[56:] if i % 2 == 0 else s for i, s in enumerate(foreground)]
(t / "fg.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(foreground)))
(t / "bg.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(5000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_best_motif_in_cluster(
    "motifs.pfm",
    "motifs.pfm",
    [{"representative": "MA0001.1 AR1", "members": ["MA0001.1 AR1", "MA0001.2 AR1b"]}],
    str(t / "fg.fa"),
    {"random": str(t / "bg.fa")},
    str(t / "genome.fa"),
    verbose=False,
)
print(result["best_motifs"])
```

- **示例输出：**

```text
['MA0001.1 AR1']
```

## 27) `dnakit.sequence.gimmemotifs_denovo_create_background()`

- **作用：** 为从头基序分析创建一个背景 FASTA 文件。
- **处理方法：** 对 `random`、`genomic`、`gc` 和 `custom`，调用 GimmeMotifs 从头预测模块的 `create_background()`。对 `promoter`，用 `PromoterFasta` 根据基因 BED 文件构建随机启动子并写出。每条前景序列对应生成 `nr_times` 条背景序列。基因组始终是 FASTA 路径。
- **API：** `gimmemotifs_denovo_create_background(background_type: 'str', fasta_file: 'str', output_path: 'str', *, genome_fasta: 'str | None' = None, annotation_bed: 'str | None' = None, size: 'int' = 200, nr_times: 'int' = 10, custom_background: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。`gc` 类型还需要 PATH 中有 `bedtools`。
- **输入：** `background_type` 为 `random`（`fasta_file` 的一阶马尔可夫模型）、`genomic`（`genome_fasta` 的随机区域）、`gc`（与 `fasta_file` GC 含量匹配的 `genome_fasta` 区域）、`promoter`（需要 `genome_fasta` 和 `annotation_bed`）或 `custom`（原样复制 `custom_background`）。`output_path` 不能已存在。`size` 和 `nr_times` 是正整数。
- **输出：** `background_type`、`sequence_count` 以及含 `background` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_denovo_create_background

out = Path(tempfile.mkdtemp()) / "background.fa"
result = gimmemotifs_denovo_create_background(
    "random", "seqs.fa", out, nr_times=3, verbose=False
)
records = out.read_text().split(">")[1:]
print(result["background_type"], result["sequence_count"], len(records))
```

- **示例输出：**

```text
random 18 18
```

## 28) `dnakit.sequence.gimmemotifs_denovo_create_backgrounds()`

- **作用：** 为一次从头预测创建预测背景和验证背景。
- **处理方法：** 把预测与验证 FASTA 文件作为 `prediction.fa` 和 `validation.fa` 复制到新的输出目录，然后调用 GimmeMotifs 的 `create_backgrounds()`。预测背景使用第一个类型（列出 `gc` 时用 `gc`）；每个类型各写出一个验证背景文件。
- **API：** `gimmemotifs_denovo_create_backgrounds(output_dir: 'str', *, prediction_fasta: 'str', validation_fasta: 'str', background_types: 'Sequence[str]' = ('random',), genome_fasta: 'str | None' = None, size: 'int' = 200, custom_background: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。`gc` 类型还需要 PATH 中有 `bedtools`。
- **输入：** `output_dir` 不能已存在。`prediction_fasta` 和 `validation_fasta` 是已存在的文件。`background_types` 是取自 `random`、`genomic`、`gc`、`custom` 的列表（`promoter` 请用 `gimmemotifs_denovo_create_background`）。`genomic` 和 `gc` 需要 `genome_fasta`；`custom` 需要 `custom_background`。
- **输出：** `background_types`、`background_files`（类型对应验证文件路径），以及 `outputs`（含 `prediction_background`，每个类型一项 `validation_<type>`）。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_denovo_create_backgrounds

t = Path(tempfile.mkdtemp())
rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


(t / "prediction.fa").write_text("".join(f">p{i}\n{seq(100)}\n" for i in range(10)))
(t / "validation.fa").write_text("".join(f">v{i}\n{seq(100)}\n" for i in range(10)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(3000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_denovo_create_backgrounds(
    t / "run",
    prediction_fasta=str(t / "prediction.fa"),
    validation_fasta=str(t / "validation.fa"),
    background_types=["random", "genomic"],
    genome_fasta=str(t / "genome.fa"),
    size=100,
    verbose=False,
)
print(sorted(result["outputs"]))
print(sorted(path.name for path in (t / "run").iterdir()))
```

- **示例输出：**

```text
['prediction_background', 'validation_genomic', 'validation_random']
['bg.genomic.fa', 'bg.random.fa', 'prediction.bg.fa', 'prediction.fa', 'validation.fa']
```

## 29) `dnakit.sequence.gimmemotifs_denovo_motif_report()`

- **作用：** 为一组从头预测基序写出文本和 HTML 报告。
- **处理方法：** 按上游报告的步骤执行：为每个基序和背景画 ROC 图，用 `MotifComparer` 找最接近的数据库基序，计算统计（`calc_stats()` 关闭 z 分数和 GC，再 `add_star()`），根据 `location_fasta` 为每个基序画定位直方图，最后写出文本报告和 HTML 报告。写出的每个文件都会返回。
- **API：** `gimmemotifs_denovo_motif_report(input_file: 'str', pfm_file: 'str', foreground_fasta: 'str', background_fastas: 'Mapping[str, str]', location_fasta: 'str', output_dir: 'str', *, genome_fasta: 'str', cutoff_fpr: 'float' = 0.9, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib。
- **输入：** `input_file` 是已存在的文件，其名称显示在报告中。`pfm_file` 含这些基序。`foreground_fasta`、`location_fasta`（以峰为中心的序列）和 `genome_fasta` 是已存在的 FASTA 文件。`background_fastas` 把背景名映射到 FASTA 路径。`output_dir` 不能已存在。`cutoff_fpr` 是直方图的得分阈值（默认 0.9）。
- **输出：** `report_folder` 和 `outputs`（相对文件路径对应路径）：`gimme.denovo.html`、`stats.<background>.txt`，以及每个基序的序列标志、`<id>_roc.<background>.png` 和 `<id>_histogram.svg`。最接近的数据库基序的序列标志也会列出。
- **示例代码：**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_denovo_motif_report

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


foreground = [seq(200) for _ in range(40)]
foreground = [s[:50] + "TGACAC" + s[56:] if i % 2 == 0 else s for i, s in enumerate(foreground)]
(t / "fg.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(foreground)))
(t / "bg.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(5000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_denovo_motif_report(
    str(t / "fg.fa"),
    "motifs.pfm",
    str(t / "fg.fa"),
    {"random": str(t / "bg.fa")},
    str(t / "fg.fa"),
    t / "report",
    genome_fasta=str(t / "genome.fa"),
    verbose=False,
)
names = sorted(Path(path).name for path in result["outputs"].values())
print("gimme.denovo.html" in names, "stats.random.txt" in names)
print([n for n in names if n.endswith("_roc.random.png")])
print([n for n in names if n.endswith("_histogram.svg")])
```

- **示例输出：**

```text
True True
['MA0001.1 AR1_roc.random.png', 'MA0001.2 AR1b_roc.random.png', 'MA0002.1 GATA_roc.random.png']
['MA0001.1 AR1_histogram.svg', 'MA0001.2 AR1b_histogram.svg', 'MA0002.1 GATA_histogram.svg']
```

## 30) `dnakit.sequence.gimmemotifs_diff_motifs()`

- **作用：** 比较若干 FASTA 集合与背景 FASTA 的基序频率，并画出差异。
- **处理方法：** 遵循 `gimme diff`：用阈值扫描每个基序，统计背景和每个前景集合中含命中的序列数；频率为 `(count + 0.01) / n`，其中 `n` 为该集合的序列数。集合名取 FASTA 文件主干。然后用这些表调用 `gimmemotifs_diff_plot()`。上游比较的是 `log2(频率 / 背景频率)` 与 `min_enrichment`，所以默认值 3 表示八倍富集，并且它会把保留的基序打印到标准输出。
- **API：** `gimmemotifs_diff_motifs(foreground_fastas: 'Sequence[str]', background_fasta: 'str', pfm_file: 'str', output_path: 'str', *, cutoff: 'float' = 0.9, min_enrichment: 'float' = 3, min_frequency: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib 和 scipy。
- **输入：** `foreground_fastas` 是非空的 FASTA 路径列表，文件主干互不相同。`background_fasta` 和 `pfm_file` 是已存在的文件。`output_path` 不能已存在。`cutoff` 是得分阈值（不小于 0）。`min_enrichment` 和 `min_frequency` 是绘图阈值。
- **输出：** `gimmemotifs_diff_plot()` 的字段（`plot_written`、`motif_ids`、`sample_names`、`frequencies`、`counts`、`background_frequencies`、`background_counts`），外加 `cutoff`、`background_size` 和 `foreground_sizes`。
- **示例代码：**

```python
import contextlib
import io
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_diff_motifs

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def plant(site, n):
    out = []
    for _ in range(n):
        s = seq(100)
        out.append(s[:40] + site + s[40 + len(site):])
    return out


sets = {"setA": plant("TGACAC", 30), "setB": plant("GCATAC", 30), "bg": [seq(100) for _ in range(60)]}
for name, seqs in sets.items():
    (t / f"{name}.fa").write_text("".join(f">{name}{i}\n{s}\n" for i, s in enumerate(seqs)))
with contextlib.redirect_stdout(io.StringIO()):
    result = gimmemotifs_diff_motifs(
        [str(t / "setA.fa"), str(t / "setB.fa")],
        str(t / "bg.fa"),
        "motifs.pfm",
        t / "diff.png",
        min_enrichment=1,
        verbose=False,
    )
print(result["plot_written"], result["background_size"], result["foreground_sizes"])
for motif, row in zip(result["motif_ids"], result["counts"]):
    print(motif, row)
```

- **示例输出：**

```text
True 60 {'setA': 30, 'setB': 30}
MA0001.1 AR1 [30.0, 4.0]
MA0002.1 GATA [6.0, 30.0]
MA0001.2 AR1b [30.0, 3.0]
```

## 31) `dnakit.sequence.gimmemotifs_diff_plot()`

- **作用：** 画出若干序列集合的基序频率与背景的对比。
- **处理方法：** 检查表的形状并调用 GimmeMotifs 的 `diff_plot()`。上游只在以下条件都满足时保留某个基序：`log2(频率 / 背景频率)` 在某个集合中超过 `min_enrichment`，某个频率超过 `min_frequency`，各集合间对数富集的差距超过 `min_difference`，且某个计数大于 2。没有基序通过时不画图（`plot_written` 为 False）；只有一个集合时差距为 0，所以在默认的 `min_difference` 下不会画图。上游也会把保留的基序打印到标准输出。
- **API：** `gimmemotifs_diff_plot(motif_ids: 'Sequence[str]', pfm_file: 'str', sample_names: 'Sequence[str]', frequencies: 'Sequence[Sequence[float]]', counts: 'Sequence[Sequence[float]]', background_frequencies: 'Sequence[float]', background_counts: 'Sequence[float]', output_path: 'str', *, min_difference: 'float' = 0, min_enrichment: 'float' = 3, min_frequency: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib。
- **输入：** `motif_ids` 是 `pfm_file` 中的 ID。`frequencies` 和 `counts` 每个基序一行、`sample_names` 的每个条目一列；`background_frequencies` 和 `background_counts` 每个基序一个值。`output_path` 不能已存在。阈值是数字（`min_frequency` 不小于 0）。
- **输出：** `plot_written`、按原样返回的输入表（`motif_ids`、`sample_names`、`frequencies`、`counts`、`background_frequencies`、`background_counts`），以及写出图时含 `plot` 的 `outputs`。
- **示例代码：**

```python
import contextlib
import io
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_diff_plot

t = Path(tempfile.mkdtemp())
with contextlib.redirect_stdout(io.StringIO()):
    result = gimmemotifs_diff_plot(
        ["MA0001.1 AR1", "MA0002.1 GATA"],
        "motifs.pfm",
        ["setA", "setB"],
        [[0.8, 0.1], [0.1, 0.7]],
        [[24, 3], [3, 21]],
        [0.1, 0.1],
        [6, 6],
        t / "diff.png",
        min_enrichment=1,
        verbose=False,
    )
    single = gimmemotifs_diff_plot(
        ["MA0001.1 AR1"], "motifs.pfm", ["setA"], [[0.8]], [[24]], [0.1], [6], t / "one.png",
        min_enrichment=1, verbose=False,
    )
print(result["plot_written"], sorted(result["outputs"]), single["plot_written"], single["outputs"])
```

- **示例输出：**

```text
True ['plot'] False {}
```

## 32) `dnakit.sequence.gimmemotifs_factor2orthogroups()`

- **作用：** 返回某个转录因子在直系同源数据库中的同源群。
- **处理方法：** 用 `strict` 策略调用 GimmeMotifs 的 `factor2orthogroups()`：在 SQLite 数据库中，把因子名或 ID（不区分大小写）与参考组装的基因匹配。medium 和 lenient 策略会到 mygene.info 查询名称，因此被拒绝。上游在进程存续期间按因子、参考组装和数据库路径缓存结果。
- **API：** `gimmemotifs_factor2orthogroups(factor: 'str', references: 'Sequence[str]', database_path: 'str', *, strategy: 'str' = 'strict', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），以及由 `gimmemotifs_load_orthogroups_in_db()` 生成的数据库。
- **输入：** `factor` 是不含引号的非空名称。`references` 是数据库中组装名的非空列表。`database_path` 是已存在的 SQLite 文件。`strategy` 必须是 `"strict"`。
- **输出：** `factor` 和 `orthogroups`（排序后的同源群 ID 列表；找不到因子或因子在上游黑名单中时为空）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_factor2orthogroups, gimmemotifs_load_orthogroups_in_db

t = Path(tempfile.mkdtemp())
(t / "results" / "Orthogroups").mkdir(parents=True)
(t / "results" / "Orthogroups" / "Orthogroups.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\n"
    "OG0000000\tGATA1|ENSG1\tgata1a|ENSDARG1, gata1b|ENSDARG2\n"
    "OG0000001\tTP53|ENSG2\ttp53|ENSDARG3\n"
)
(t / "results" / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\nOG0000002\t.|ENSG9\t\n"
)
database = t / "orthologs.sqlite"
gimmemotifs_load_orthogroups_in_db(database, ["human", "fish"], str(t / "results"), verbose=False)
for factor in ("gata1", "ENSG2", "unknown"):
    result = gimmemotifs_factor2orthogroups(factor, ["human"], str(database), verbose=False)
    print(result["factor"], result["orthogroups"])
```

- **示例输出：**

```text
gata1 ['OG0000000']
ENSG2 ['OG0000001']
unknown []
```

## 33) `dnakit.sequence.gimmemotifs_filter_significant_motifs()`

- **作用：** 保留统计量达到最低值的从头预测基序。
- **处理方法：** 读取基序，把每个基序的统计量以 `<id>_<共有序列>` 为键（GimmeMotifs 使用的键）组织起来，并针对所选背景类型调用 GimmeMotifs 的 `filter_significant_motifs()`。显著的基序写成 PFM 文本。只有 `metrics` 中每个统计量都达到各自最低值的基序才会保留。
- **API：** `gimmemotifs_filter_significant_motifs(motif_file: 'str', output_path: 'str', *, background_type: 'str', stats: 'Mapping[str, Mapping[str, Mapping[str, float]]]', metrics: 'Mapping[str, float] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`）。
- **输入：** `motif_file` 是已存在的 PFM 文件。`output_path` 不能已存在。`background_type` 是非空字符串，例如 `random`。`stats` 把基序 ID 映射到背景类型再映射到统计量（`max_enrichment`、`roc_auc`、`enr_at_fpr` 等）；没有对应条目的基序会被丢弃。`metrics` 把统计量名映射到最低值；为 None 时使用上游默认值（`max_enrichment` 3、`roc_auc` 0.55、`enr_at_fpr` 0.55）。
- **输出：** `significant_motifs`（保留的 ID）、`background_type` 以及含 `motifs` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_filter_significant_motifs

out = Path(tempfile.mkdtemp()) / "significant.pfm"
stats = {
    "MA0001.1 AR1": {"random": {"max_enrichment": 5.0, "roc_auc": 0.8, "enr_at_fpr": 0.9}},
    "MA0002.1 GATA": {"random": {"max_enrichment": 1.0, "roc_auc": 0.5, "enr_at_fpr": 0.2}},
}
result = gimmemotifs_filter_significant_motifs(
    "motifs.pfm", out, background_type="random", stats=stats, verbose=False
)
print(result["significant_motifs"], out.read_text().count(">"))
```

- **示例输出：**

```text
['MA0001.1 AR1'] 1
```

## 34) `dnakit.sequence.gimmemotifs_load_orthogroups_in_db()`

- **作用：** 把 OrthoFinder 的同源群存入新的 SQLite 数据库。
- **处理方法：** 检查两个 OrthoFinder 表存在且每个组装都有 `<assembly>.pep` 列，然后调用 GimmeMotifs 的 `load_orthogroups_in_db()`。它创建 `orthogroups`、`assemblies` 和 `genes` 三张表，并为基因名与 ID 建索引。基因以 `gene_name|gene_id` 的形式写入（名称缺失时为 `.`）。随后函数统计行数。
- **API：** `gimmemotifs_load_orthogroups_in_db(database_path: 'str', assemblies: 'Sequence[str]', orthofinder_results_dir: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），以及 pandas。不会运行 OrthoFinder 本身。
- **输入：** `database_path` 不能已存在（上游会替换它）。`assemblies` 是非空的名称列表。`orthofinder_results_dir` 含 `Orthogroups/Orthogroups.tsv` 和 `Orthogroups/Orthogroups_UnassignedGenes.tsv`。
- **输出：** `orthogroup_count`（含 `UNASSIGNED` 组）、每个组装的 `gene_counts`，以及含 `database` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_load_orthogroups_in_db

t = Path(tempfile.mkdtemp())
(t / "results" / "Orthogroups").mkdir(parents=True)
(t / "results" / "Orthogroups" / "Orthogroups.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\n"
    "OG0000000\tGATA1|ENSG1\tgata1a|ENSDARG1, gata1b|ENSDARG2\n"
    "OG0000001\tTP53|ENSG2\ttp53|ENSDARG3\n"
)
(t / "results" / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\nOG0000002\t.|ENSG9\t\n"
)
result = gimmemotifs_load_orthogroups_in_db(
    t / "orthologs.sqlite", ["human", "fish"], str(t / "results"), verbose=False
)
print(result["orthogroup_count"], result["gene_counts"], sorted(result["outputs"]))
```

- **示例输出：**

```text
3 {'human': 3, 'fish': 3} ['database']
```

## 35) `dnakit.sequence.gimmemotifs_make_motif2factors()`

- **作用：** 根据直系同源数据库，为新的基因组组装写出基序到因子的对应表。
- **处理方法：** 对每个基序，在 `database_references` 组装中找到其因子所在的同源群，并把 `new_reference` 中属于这些同源群的每个基因都写成该基序的因子（GimmeMotifs 的 `make_motif2factors()`）。没有任何同源基因的基序会得到 `NO ORTHOLOGS FOUND` 一行。基序的 PFM 文件写在对应表旁边。因为其他策略要联网，只允许 `strict` 策略。
- **API：** `gimmemotifs_make_motif2factors(output_prefix: 'str', *, motif_file: 'str', motif_factors: 'Mapping[str, Sequence[str]]', database_path: 'str', new_reference: 'str', database_references: 'Sequence[str]', strategy: 'str' = 'strict', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），以及由 `gimmemotifs_load_orthogroups_in_db()` 生成的同源数据库。
- **输入：** `output_prefix` 是路径前缀；`<prefix>.motif2factors.txt` 和 `<prefix>.pfm` 不能已存在。`motif_file` 是 PFM 文件。`motif_factors` 把基序 ID 映射到因子名列表。`database_path` 是已存在的数据库；`new_reference` 和 `database_references` 是其中的组装名。`strategy` 必须是 `"strict"`。
- **输出：** `motif_count` 以及含 `motif2factors` 和 `pfm` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_load_orthogroups_in_db, gimmemotifs_make_motif2factors

t = Path(tempfile.mkdtemp())
(t / "results" / "Orthogroups").mkdir(parents=True)
(t / "results" / "Orthogroups" / "Orthogroups.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\n"
    "OG0000000\tGATA1|ENSG1\tgata1a|ENSDARG1, gata1b|ENSDARG2\n"
    "OG0000001\tTP53|ENSG2\ttp53|ENSDARG3\n"
)
(t / "results" / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(
    "Orthogroup\thuman.pep\tfish.pep\nOG0000002\t.|ENSG9\t\n"
)
database = t / "orthologs.sqlite"
gimmemotifs_load_orthogroups_in_db(database, ["human", "fish"], str(t / "results"), verbose=False)
result = gimmemotifs_make_motif2factors(
    str(t / "fish"),
    motif_file="motifs.pfm",
    motif_factors={"MA0001.1 AR1": ["GATA1"], "MA0002.1 GATA": ["TP53"], "MA0001.2 AR1b": ["XYZ"]},
    database_path=str(database),
    new_reference="fish",
    database_references=["human"],
    verbose=False,
)
print(result["motif_count"], sorted(result["outputs"]))
for line in sorted((t / "fish.motif2factors.txt").read_text().splitlines()[1:]):
    print(line.replace("\t", " | "))
```

- **示例输出：**

```text
3 ['motif2factors', 'pfm']
MA0001.1 AR1 | gata1a | Orthologs | N
MA0001.1 AR1 | gata1b | Orthologs | N
MA0001.2 AR1b | NO ORTHOLOGS FOUND | Orthologs | N
MA0002.1 GATA | tp53 | Orthologs | N
```

## 36) `dnakit.sequence.gimmemotifs_motif2factor_from_orthologs()`

- **作用：** 根据本地基因组和注释，通过直系同源推断为新的组装构建基序到因子的对应表。
- **处理方法：** 检查 PATH 中有 `gffread` 和 `orthofinder`；然后对每个组装，把基因组和 GTF 复制到 `intermediate/annotated/<name>/` 并运行 GimmeMotifs 的 `annot2primpep()`；在肽段目录上运行 OrthoFinder；把同源群载入 SQLite 数据库；并对 `new_references` 中的每个名称用 `strict` 策略调用 `make_motif2factors()`。不下载任何内容。
- **API：** `gimmemotifs_motif2factor_from_orthologs(genomes: 'Mapping[str, Mapping[str, str]]', database_references: 'Sequence[str]', new_references: 'Sequence[str]', motif_file: 'str', motif_factors: 'Mapping[str, Sequence[str]]', output_dir: 'str', *, threads: 'int' = 1, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），并且 PATH 中有外部程序 `gffread` 和 `orthofinder`（例如 `conda install -c bioconda gffread orthofinder`）。下面的示例在缺少它们的情况下运行并展示报错；验证本页所用的环境没有安装这两个程序。
- **输入：** `genomes` 把每个组装名映射到 `{"fasta": 路径, "annotation": 路径}`。`database_references` 和 `new_references` 是名称列表，名称必须是 `genomes` 的键。`motif_file` 是 PFM 文件，`motif_factors` 把基序 ID 映射到因子名。`output_dir` 不能已存在。`threads` 是正整数。
- **输出：** `orthofinder_results`（OrthoFinder 结果目录）以及 `outputs`，其中每个新参考对应一项 `<reference>_motif2factors`，路径为 `output_dir/<reference>.<基序文件名>.motif2factors.txt`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence import gimmemotifs_motif2factor_from_orthologs

out = Path(tempfile.mkdtemp()) / "motif2factor"
genomes = {name: {"fasta": "seqs.fa", "annotation": "motifs.pfm"} for name in ("human", "fish")}
try:
    gimmemotifs_motif2factor_from_orthologs(
        genomes, ["human"], ["fish"], "motifs.pfm", {"MA0001.1 AR1": ["GATA1"]}, out, verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(" Hint")[0])
print(out.exists())
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] The gffread program was not found on PATH.
False
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_motif2factor_from_orthologs

result = gimmemotifs_motif2factor_from_orthologs(
    {
        "human": {"fasta": "human.fa", "annotation": "human.gtf"},
        "mouse": {"fasta": "mouse.fa", "annotation": "mouse.gtf"},
    },
    database_references=["human"],
    new_references=["mouse"],
    motif_file="motifs.pfm",
    motif_factors={"MA0001.1 AR1": ["GATA1"]},
    output_dir="motif2factor_out",
    threads=4,
)
print(result["outputs"])
```

## 37) `dnakit.sequence.gimmemotifs_motif_tree_plot()`

- **作用：** 画出基序聚类树，分支颜色由每个基序的数值决定。
- **处理方法：** 检查 ete3 能否导入（在写出任何内容之前），验证 Newick 文本和每个基序的数值，构建每个基序一列的 pandas 表，并调用 GimmeMotifs 的 `motif_tree_plot()`。分支颜色取其下各叶子数值的平均值。
- **API：** `gimmemotifs_motif_tree_plot(output_path: 'str', newick_tree: 'str', leaf_values: 'Mapping[str, Sequence[float]]', *, circle: 'bool' = True, vmin: 'float | None' = None, vmax: 'float | None' = None, dpi: 'int' = 300, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs` 和 `ete3`（都在 extra `dnakit[bio-gimmemotifs]` 中）。验证本页所用的环境没有安装 ete3，所以示例展示在写出任何文件之前就抛出的错误。
- **输入：** `output_path` 不能已存在。`newick_tree` 是以 `;` 结尾的 Newick 文本，叶子为基序 ID。`leaf_values` 把每个基序 ID 映射到数字列表，所有列表长度相同。`circle` 选择环形布局。`vmin` 和 `vmax` 设置颜色范围。`dpi` 是正整数。
- **输出：** `leaf_count` 以及含 `plot` 的 `outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence import gimmemotifs_motif_tree_plot

out = Path(tempfile.mkdtemp()) / "tree.png"
try:
    gimmemotifs_motif_tree_plot(
        out, "((m1,m2),m3);", {"m1": [1.0, 2.0], "m2": [1.5, 2.5], "m3": [-1.0, 0.0]}, verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(" Hint")[0])
print(out.exists())
```

- **示例输出：**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] ete3 could not be imported: No module named 'ete3'
False
```

安装依赖后的用法：

```python
from dnakit.sequence import gimmemotifs_motif_tree_plot

result = gimmemotifs_motif_tree_plot(
    "tree.png",
    "((m1,m2),m3);",
    {"m1": [1.0, 2.0], "m2": [1.5, 2.5], "m3": [-1.0, 0.0]},
    circle=True,
    dpi=150,
)
print(result["leaf_count"], result["outputs"]["plot"])
```

## 38) `dnakit.sequence.gimmemotifs_roc_html_report()`

- **作用：** 为基序富集表写出 HTML 报告。
- **处理方法：** 读取该表，加入 Benjamini-Hochberg 校正后的 P 值及其负 log10，保留校正后 P 值不高于 `threshold` 的基序，把它们的序列标志画到 `<output_dir>/logos`，并写出 HTML 文件（GimmeMotifs 的 `roc_html_report()`）。
- **API：** `gimmemotifs_roc_html_report(output_dir: 'str', input_table: 'str', pfm_file: 'str', output_name: 'str', *, threshold: 'float' = 0.01, use_motifs: 'Sequence[str] | None' = None, link_matches: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），以及 pandas。
- **输入：** `output_dir` 不能已存在。`input_table` 是已完成运行的制表符分隔富集表（第一列为基序 ID，其余列为 `# matches`、`% matches input`、`# matches background`、`% matches background`、`P-value`、`log10 P-value`、`ROC AUC`、`PR AUC`、`Enr. at 1% FPR`、`Recall at 10% FDR`）。`pfm_file` 含这些基序。`output_name` 是普通文件名。`threshold` 是 [0, 1] 内的数。`use_motifs` 限定 ID。`link_matches` 把匹配数链接到扫描得到的 BED 文件。
- **输出：** `threshold` 以及 `outputs`：`report`，加上目录下的每个文件（HTML 文件和序列标志，按相对路径）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_roc_html_report

t = Path(tempfile.mkdtemp())
header = (
    "Motif\t# matches\t% matches input\t# matches background\t% matches background\t"
    "P-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n"
)
(t / "roc.txt").write_text(
    header
    + "MA0001.1 AR1\t30\t75.0\t8\t10.0\t1e-12\t12.0\t0.90\t0.85\t6.0\t0.50\n"
    + "MA0002.1 GATA\t20\t50.0\t10\t12.5\t1e-04\t4.0\t0.70\t0.65\t3.0\t0.20\n"
    + "MA0001.2 AR1b\t28\t70.0\t8\t10.0\t1e-10\t10.0\t0.88\t0.80\t5.0\t0.45\n"
)
result = gimmemotifs_roc_html_report(
    t / "report", str(t / "roc.txt"), "motifs.pfm", "report.html", threshold=1e-6, verbose=False
)
print(sorted(result["outputs"]))
html = (t / "report" / "report.html").read_text()
print("MA0001.1 AR1" in html, "MA0002.1 GATA" in html)
```

- **示例输出：**

```text
['logos/MA0001_1_AR1.png', 'logos/MA0001_2_AR1b.png', 'report', 'report.html']
True False
```

## 39) `dnakit.sequence.gimmemotifs_roc_plots()`

- **作用：** 为每个基序和每个背景集合画一条 ROC 曲线。
- **处理方法：** 对每个基序和每个背景 FASTA，在前景和背景上计算该基序的 ROC 值（`calc_stats()` 关闭 z 分数和 GC 选项，因为默认选项会用 `str` 对基序做哈希，而 xxhash 3 及更高版本拒绝 `str`），并用 GimmeMotifs 的报告函数画出曲线。函数会检查每一对都写出了一个文件。
- **API：** `gimmemotifs_roc_plots(pfm_file: 'str', foreground_fasta: 'str', background_fastas: 'Mapping[str, str]', output_dir: 'str', *, genome_fasta: 'str', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 matplotlib。
- **输入：** `pfm_file`、`foreground_fasta` 和 `genome_fasta` 是已存在的文件。`background_fastas` 把背景名映射到 FASTA 路径。`output_dir` 不能已存在。
- **输出：** `motif_count`、`background_count`，以及 `outputs`（`<motif id>_roc.<background>.png` 对应其路径；文件位于 `output_dir/images`）。
- **示例代码：**

```python
import os
import random
import tempfile
from pathlib import Path

t = Path(tempfile.mkdtemp())
os.environ["XDG_CACHE_HOME"] = str(t / "cache")

from dnakit.sequence import gimmemotifs_roc_plots

rng = random.Random(1)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


foreground = [seq(200) for _ in range(40)]
foreground = [s[:50] + "TGACAC" + s[56:] if i % 2 == 0 else s for i, s in enumerate(foreground)]
(t / "fg.fa").write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(foreground)))
(t / "bg.fa").write_text("".join(f">b{i}\n{seq(200)}\n" for i in range(80)))
(t / "genome.fa").write_text("".join(f">{c}\n{seq(5000)}\n" for c in ("chr1", "chr2")))
result = gimmemotifs_roc_plots(
    "motifs.pfm",
    str(t / "fg.fa"),
    {"random": str(t / "bg.fa")},
    t / "roc",
    genome_fasta=str(t / "genome.fa"),
    verbose=False,
)
print(result["motif_count"], result["background_count"])
print(sorted(result["outputs"]))
```

- **示例输出：**

```text
3 1
['MA0001.1 AR1_roc.random.png', 'MA0001.2 AR1b_roc.random.png', 'MA0002.1 GATA_roc.random.png']
```

## 40) `dnakit.sequence.gimmemotifs_select_nonredundant_motifs()`

- **作用：** 根据 ROC 报告和两个基序得分表选出一组非冗余基序。
- **处理方法：** 调用 GimmeMotifs 的 `select_nonredundant_motifs()`。`Enr. at 1% FPR` 低于 2 的基序被丢弃，其余按报告中四列的平均秩排序，与更好基序相似（seqcor 不低于 0.7）的基序被去掉；如果剩下至少两个基序，再用 L1 逻辑回归做递归特征消除，保留能提高平均精确率的基序。
- **API：** `gimmemotifs_select_nonredundant_motifs(roc_report: 'str', pfm_file: 'str', foreground_table: 'str', background_table: 'str', *, tolerance: 'float' = 0.001, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `gimmemotifs`（extra `dnakit[bio-gimmemotifs]`），它使用 scikit-learn。
- **输入：** `roc_report` 是制表符分隔的表，第一列为基序 ID，另有 `ROC AUC`、`PR AUC`、`Enr. at 1% FPR` 和 `Recall at 10% FDR` 列。`pfm_file` 含这些基序。`foreground_table` 和 `background_table` 是制表符分隔的得分表，每条序列一行、每个基序 ID 一列。`tolerance` 不小于 0。
- **输出：** `selected_motifs`，即选出的基序 ID。
- **示例代码：**

```python
import random
import tempfile
from pathlib import Path

from dnakit.sequence import gimmemotifs_select_nonredundant_motifs

t = Path(tempfile.mkdtemp())
header = (
    "Motif\t# matches\t% matches input\t# matches background\t% matches background\t"
    "P-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n"
)
(t / "roc.txt").write_text(
    header
    + "MA0001.1 AR1\t30\t75.0\t8\t10.0\t1e-12\t12.0\t0.90\t0.85\t6.0\t0.50\n"
    + "MA0002.1 GATA\t20\t50.0\t10\t12.5\t1e-04\t4.0\t0.70\t0.65\t3.0\t0.20\n"
    + "MA0001.2 AR1b\t28\t70.0\t8\t10.0\t1e-10\t10.0\t0.88\t0.80\t5.0\t0.45\n"
)
rng = random.Random(3)


def table(n, shift):
    rows = ["\tMA0001.1 AR1\tMA0002.1 GATA\tMA0001.2 AR1b"]
    for i in range(n):
        a = rng.gauss(shift, 1)
        rows.append(f"s{i}\t{a:.3f}\t{rng.gauss(shift / 2, 1):.3f}\t{a + rng.gauss(0, 0.1):.3f}")
    return "\n".join(rows) + "\n"


(t / "fg.tsv").write_text(table(60, 2))
(t / "bg.tsv").write_text(table(60, 0))
result = gimmemotifs_select_nonredundant_motifs(
    str(t / "roc.txt"), "motifs.pfm", str(t / "fg.tsv"), str(t / "bg.tsv"), verbose=False
)
print(result["selected_motifs"])
```

- **示例输出：**

```text
['MA0001.1 AR1']
```
