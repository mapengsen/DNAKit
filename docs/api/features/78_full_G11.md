# 表观组函数：methylpy、MACS3 与 deepTools

这一组函数位于 `dnakit.epigenome`，是 methylpy、MACS3 与 deepTools 的 Python 接口的封装。每个函数直接调用上游的 Python API，不运行命令行，也不联网下载数据；可选依赖只在调用时导入，缺少时报错并给出安装命令。methylpy 与 deepTools 在输入不合法时会调用 `sys.exit`，这些退出被转换为 `ConfigurationError`。

本批共封装 32 个函数：methylpy 10 个、MACS3 12 个、deepTools 10 个。示例从模块文件导入（如 `dnakit.epigenome`），导出到 `dnakit.epigenome` 后同样可以用 `dnakit.epigenome.<名称>()` 调用。未封装的函数（命令行入口、绘图、依赖 samtools、bowtie2、cutadapt、bgzip、tabix 或 wigToBigWig 的流程、需要全局状态的 GC 偏差函数）见账本，原因逐条记录。

## 示例数据

示例数据是小型文本文件，下面的示例按顺序在同一目录中运行，后面的示例会使用前面生成的文件。

带索引的 FASTA（读取染色体序列）。

**示例文件：** `genome.fa`

```text
>chr1
ACGTACGTNNACGT
>chr2
GGCC
```

FASTQ（C 到 T 转换）。

**示例文件：** `reads.fq`

```text
@r1
ACGTCC
+
IIIIII
@r2
TTCG
+
IIII
```

30 bp 参考序列（mCH 过滤）。

**示例文件：** `ref.fa`

```text
>chr1
ACAACAACAACAACAACAACAAAAAAAAAA
```

两条比对 read，一条高甲基化、一条未甲基化（SAM 文本，示例中转换为 BAM）。

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:30
r1	0	chr1	1	60	20M	*	0	0	ACAACAACAACAACAACAAC	IIIIIIIIIIIIIIIIIIII
r2	0	chr1	1	60	20M	*	0	0	ATAATAATAATAATAATAAT	IIIIIIIIIIIIIIIIIIII
```

按第 7 列 P 值升序的二项检验结果。

**示例文件：** `pvalues.tsv`

```text
chr1	10	+	CGN	8	10	0.0001
chr1	40	+	CHH	5	9	0.001
chr1	20	+	CGN	1	10	0.2
chr1	30	+	CHH	0	9	0.6
```

按染色体连续排列的 allc 表。

**示例文件：** `allc.tsv`

```text
chr1	10	+	CGN	8	10	1
chr1	20	+	CGN	1	10	1
chr2	30	+	CHH	0	9	1
```

两条 BED 格式 tag。

**示例文件：** `tags.bed`

```text
chr1	10	30	r1	0	+
chr1	20	40	r2	0	-
```

两段 bedGraph 信号。

**示例文件：** `signal.bedgraph`

```text
chr1	0	10	1.5
chr1	10	20	2
```

**示例文件（genome.2bit）**：由下面的示例用 base64 写入，序列为 `ACGTGGCCAANNGGTTACGT`。

## methylpy 函数（bisulfite 甲基化）

### 1) `dnakit.epigenome.expand_methylation_context()`

- **作用：** 把带 IUPAC 符号的胞嘧啶上下文（如 `CGN`）展开成 allc 表中可能出现的具体上下文字符串。
- **处理方法：** 调用 methylpy 的 `expand_nucleotide_code`，结果排序后返回。methylpy 也会保留模式中的单字母代码，本函数同样保留。
- **API：** `dnakit.epigenome.expand_methylation_context(mc_type: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `mc_type`：1 到 6 个 IUPAC 核苷酸代码，例如 `"CGN"`、`"CHH"`。
- **输出：** `contexts` 为排序后的字符串列表，`count` 为其长度。
- **示例代码：**

```python
from dnakit.epigenome import expand_methylation_context

result = expand_methylation_context("CGN", verbose=False)
print(result["count"])
print("CGN" in result["contexts"], "CHH" in result["contexts"])
```

- **示例输出：**

```text
38
True True
```

### 2) `dnakit.epigenome.encode_read_c_positions()`

- **作用：** 记录转换后的 read 中胞嘧啶的位置，用可打印字符编码，供后续解码使用。
- **处理方法：** 调用 methylpy 的 `encode_c_positions`。每个字符保存到下一个 C 的间隔，所以编码长度等于 C 的个数。
- **API：** `dnakit.epigenome.encode_read_c_positions(seq: 'str', *, is_read2: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `seq`：A、C、G、T、N 组成的序列；`is_read2`：是否为第二端 read。
- **输出：** `encoding` 为编码字符串，其余字段回显输入。
- **示例代码：**

```python
from dnakit.epigenome import encode_read_c_positions

encoded = encode_read_c_positions("ACGTCCGA", verbose=False)
print(encoded["encoding"])
print(len(encoded["encoding"]))
```

- **示例输出：**

```text
#$"
3
```

### 3) `dnakit.epigenome.decode_read_c_positions()`

- **作用：** 用 `encode_read_c_positions` 的编码把转换后的 read 恢复为原始胞嘧啶。
- **处理方法：** 调用 methylpy 的 `decode_c_positions`。正链恢复为 C，负链恢复为 G，因为 SAM 中序列总是按正链表示。
- **API：** `dnakit.epigenome.decode_read_c_positions(seq: 'str', indexes: 'str', strand: 'str', *, is_read2: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `seq`：转换后的序列；`indexes`：编码字符串（可为空）；`strand`：`"+"` 或 `"-"`。
- **输出：** `sequence` 为恢复后的序列。
- **示例代码：**

```python
from dnakit.epigenome import decode_read_c_positions, encode_read_c_positions

encoding = encode_read_c_positions("ACGTCCGA", verbose=False)["encoding"]
decoded = decode_read_c_positions("ATGTTTGA", encoding, "+", verbose=False)
print(decoded["sequence"])
```

- **示例输出：**

```text
ACGTCCGA
```

### 4) `dnakit.epigenome.count_bisulfite_basecalls()`

- **作用：** 从 pileup 碱基字符统计某个参考碱基的一致调用、不一致调用，以及 C 的转换情况。
- **处理方法：** 调用 methylpy 的 `analyze_read_basecalls`。`.` 和 `,` 表示与参考一致，`T`（或 `a`）表示 C 被转换；插入缺失被跳过。
- **API：** `dnakit.epigenome.count_bisulfite_basecalls(ref_base: 'str', read_bases: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `ref_base`：`A`、`C`、`G`、`T` 或 `N`；`read_bases`：pileup 碱基字符串。
- **输出：** `consistent_calls`、`inconsistent_calls`；参考为 C 或 G 时还有 `unconverted_c` 和 `converted_c`，否则为 `None`。
- **示例代码：**

```python
from dnakit.epigenome import count_bisulfite_basecalls

result = count_bisulfite_basecalls("C", ".,,TTA", verbose=False)
print(result["consistent_calls"], result["inconsistent_calls"])
print(result["unconverted_c"], result["converted_c"])
```

- **示例输出：**

```text
3 1
1 2
```

### 5) `dnakit.epigenome.read_fasta_chromosome()`

- **作用：** 从带索引的 FASTA 文件读取一条染色体的完整序列。
- **处理方法：** 先读取 `.fai` 索引检查长度，再调用 methylpy 的 `get_chromosome_sequence`；读到的长度与索引不一致时报错。
- **API：** `dnakit.epigenome.read_fasta_chromosome(fasta_path: 'str | Path', chrom: 'str', *, max_length: 'int' = 50000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** FASTA 文件及其同名 `.fai` 索引（可用 `pysam.faidx` 或 `samtools faidx` 生成）；`chrom`；`max_length` 为长度上限。
- **输出：** `sequence`、`length`，以及 `provenance`。
- **示例代码：**

```python
import pysam
from dnakit.epigenome import read_fasta_chromosome

pysam.faidx("genome.fa")
result = read_fasta_chromosome("genome.fa", "chr1", verbose=False)
print(result["length"], result["sequence"])
```

- **示例输出：**

```text
14 ACGTACGTNNACGT
```

### 6) `dnakit.epigenome.convert_fastq_to_bisulfite()`

- **作用：** 把 FASTQ 中每个 read 的 C 转换为 T，用于比对到亚硫酸氢盐转换后的参考序列。
- **处理方法：** 调用 methylpy 的 `convert_reads`。原始 C 的位置写入 read 名称，之后可用编码函数恢复。输出文件必须不存在。
- **API：** `dnakit.epigenome.convert_fastq_to_bisulfite(fastq_path: 'str | Path', output_path: 'str | Path', *, buffer_line_number: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `fastq_path`：FASTQ 文件；`output_path`：新的输出文件；`buffer_line_number`：缓冲行数。
- **输出：** `reads` 为 read 数；输出文件保存转换后的 FASTQ。
- **示例代码：**

```python
from pathlib import Path
from dnakit.epigenome import convert_fastq_to_bisulfite

Path("converted.fq").unlink(missing_ok=True)
result = convert_fastq_to_bisulfite("reads.fq", "converted.fq", verbose=False)
print(result["reads"])
print(Path("converted.fq").read_text().splitlines()[1])
```

- **示例输出：**

```text
2
ATGTTT
```

### 7) `dnakit.epigenome.filter_bam_by_mch_level()`

- **作用：** 删除非 CpG 甲基化水平（mCH）过高的比对 read。
- **处理方法：** 调用 methylpy 的 `bam_quality_mch_filter`。一条 read 的 CH 位点数至少为 `min_ch` 且 mCH/(mCH+uCH) 不低于 `max_mch_level` 时被删除。methylpy 接收 `min_mapq` 但不使用，本函数不暴露该参数。
- **API：** `dnakit.epigenome.filter_bam_by_mch_level(bam_path: 'str | Path', output_bam: 'str | Path', reference_fasta: 'str | Path', *, min_ch: 'int' = 3, max_mch_level: 'float' = 0.7, buffer_line_number: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 已排序并建索引的 BAM（或可直接读取的比对文件）、带 `.fai` 索引的参考 FASTA、新的 `.bam` 输出路径。
- **输出：** `input_reads`、`kept_reads`、`removed_reads`。
- **示例代码：**

```python
import pysam
from pathlib import Path
for name in ("reads.unsorted.bam", "reads.bam", "reads.bam.bai"):
    Path(name).unlink(missing_ok=True)
pysam.view("-b", "-o", "reads.unsorted.bam", "reads.sam", catch_stdout=False)
pysam.sort("-o", "reads.bam", "reads.unsorted.bam")
pysam.index("reads.bam")
import pysam
from pathlib import Path
from dnakit.epigenome import filter_bam_by_mch_level

pysam.faidx("ref.fa")
Path("filtered.bam").unlink(missing_ok=True)
result = filter_bam_by_mch_level("reads.bam", "filtered.bam", "ref.fa", verbose=False)
print(result["input_reads"], result["kept_reads"], result["removed_reads"])
```

- **示例输出：**

```text
2 1 1
```

### 8) `dnakit.epigenome.benjamini_hochberg_methylation_cutoffs()`

- **作用：** 在 Benjamini-Hochberg 错误发现率（FDR）下，求每个胞嘧啶上下文的 P 值阈值。
- **处理方法：** 调用 methylpy 的 `benjamini_hochberg_correction_call_methylated_sites`。输入必须按 P 值升序排列（第 7 列），函数先检查排序和每个上下文的检验次数。
- **API：** `dnakit.epigenome.benjamini_hochberg_methylation_cutoffs(binomial_paths: 'Sequence[str | Path]', mc_class_counts: 'Mapping[str, int]', *, sig_cutoff: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 二项检验结果文件列表；`mc_class_counts`：每个上下文的检验数；`sig_cutoff`：FDR 阈值。
- **输出：** `cutoffs` 为每个上下文的 P 值阈值（无显著位点时为 0），`tested_sites` 为读取的位点数。
- **示例代码：**

```python
from dnakit.epigenome import benjamini_hochberg_methylation_cutoffs

result = benjamini_hochberg_methylation_cutoffs(
    ["pvalues.tsv"], {"CGN": 2, "CHH": 2}, sig_cutoff=0.5, verbose=False
)
print(result["tested_sites"])
print({key: round(value, 6) for key, value in result["cutoffs"].items()})
```

- **示例输出：**

```text
4
{'CGN': 0.2, 'CHH': 0.6}
```

### 9) `dnakit.epigenome.index_allc_table()`

- **作用：** 为 allc 表建立染色体起始偏移索引，便于按染色体跳转。
- **处理方法：** 调用 methylpy 的 `index_allc_file`，索引文件 `<allc>.idx` 写在输入文件旁边，然后读取每条染色体的起始字节偏移。要求同一染色体的记录连续。
- **API：** `dnakit.epigenome.index_allc_table(allc_path: 'str | Path', *, reindex: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 按染色体排序的 allc 文件；`reindex` 为真时强制重建索引。
- **输出：** `index` 为索引路径，`chromosomes` 为每条染色体的 `chrom` 与 `offset`。
- **示例代码：**

```python
from pathlib import Path
from dnakit.epigenome import index_allc_table

Path("allc.tsv.idx").unlink(missing_ok=True)
result = index_allc_table("allc.tsv", reindex=True, verbose=False)
for item in result["chromosomes"]:
    print(item["chrom"], item["offset"])
```

- **示例输出：**

```text
chr1 0
chr2 42
```

### 10) `dnakit.epigenome.read_allc_index()`

- **作用：** 读取 allc 表索引中记录的染色体偏移。
- **处理方法：** 读取 `<allc>.idx` 并调用 methylpy 的 `read_allc_index`。索引不存在时提示先运行 `index_allc_table()`。
- **API：** `dnakit.epigenome.read_allc_index(allc_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** allc 文件路径（其 `.idx` 索引必须已存在）。
- **输出：** `chromosomes` 为按染色体名排序的 `chrom` 与 `offset` 列表。
- **示例代码：**

```python
from dnakit.epigenome import read_allc_index

result = read_allc_index("allc.tsv", verbose=False)
print([(item["chrom"], item["offset"]) for item in result["chromosomes"]])
```

- **示例输出：**

```text
[('chr1', 0), ('chr2', 42)]
```

## MACS3 函数（信号、峰与变异统计）

### 11) `dnakit.epigenome.build_paired_fragment_pileup()`

- **作用：** 由双端片段的左右坐标计算覆盖度（pileup），用于 ATAC-seq 等片段级信号。
- **处理方法：** 调用 MACS3 的 `pileup_from_LR_as_list`。输出是 `[位置, 数值]` 的分段覆盖：每个位置是一段的终点，数值是该段的覆盖。
- **API：** `dnakit.epigenome.build_paired_fragment_pileup(left_positions: 'Sequence[int]', right_positions: 'Sequence[int]', *, scale_factor: 'float' = 1.0, baseline_value: 'float' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `left_positions` 与 `right_positions`：等长的整数列表，每个片段要求右端大于左端。
- **输出：** `positions` 与 `values`（与输入等长）、`fragments`。
- **示例代码：**

```python
from dnakit.epigenome import build_paired_fragment_pileup

result = build_paired_fragment_pileup([10, 20, 25], [30, 40, 35], verbose=False)
print(result["positions"])
print(result["values"])
```

- **示例输出：**

```text
[10, 20, 25, 30, 35, 40]
[0.0, 1.0, 2.0, 3.0, 2.0, 1.0]
```

### 12) `dnakit.epigenome.build_shifted_read_pileup()`

- **作用：** 把单端 read 按链方向和偏移扩展为片段，再计算覆盖度（例如 ATAC-seq 的 Tn5 偏移）。
- **处理方法：** 调用 MACS3 的 `pileup_from_PN_shifted`。正链 read 覆盖 `[p - five_shift, p + three_shift)`，负链 read 覆盖 `[p - three_shift, p + five_shift)`，坐标截断在染色体长度内。
- **API：** `dnakit.epigenome.build_shifted_read_pileup(plus_positions: 'Sequence[int]', minus_positions: 'Sequence[int]', *, five_shift: 'int', three_shift: 'int', chrom_length: 'int', scale_factor: 'float' = 1.0, baseline_value: 'float' = 0.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 正链与负链 5' 端位置列表；`five_shift`、`three_shift`；`chrom_length`。
- **输出：** `positions`、`values`，以及 read 数和偏移。
- **示例代码：**

```python
from dnakit.epigenome import build_shifted_read_pileup

result = build_shifted_read_pileup(
    [10, 12, 30], [40, 42], five_shift=0, three_shift=10, chrom_length=100, verbose=False
)
print(result["positions"])
print(result["values"])
```

- **示例输出：**

```text
[10, 12, 20, 22, 30, 32, 40, 42]
[0.0, 1.0, 2.0, 1.0, 0.0, 2.0, 3.0, 1.0]
```

### 13) `dnakit.epigenome.find_signal_maxima()`

- **作用：** 在一维信号（如覆盖度曲线）中找局部极大值，即候选峰顶。
- **处理方法：** 调用 MACS3 的 `SignalProcessing.maxima`，先做二阶 Savitzky-Golay 平滑，再取局部极大值。`window_size` 必须为奇数且不超过信号长度。
- **API：** `dnakit.epigenome.find_signal_maxima(signal: 'Sequence[float]', *, window_size: 'int' = 51, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 数值列表（至少 3 个值）与奇数窗口 `window_size`。
- **输出：** `maxima` 为从 0 开始的索引列表。
- **示例代码：**

```python
from dnakit.epigenome import find_signal_maxima

signal = [0, 1, 2, 3, 2, 1, 0, 1, 4, 6, 4, 1, 0, 0, 0]
print(find_signal_maxima(signal, window_size=3, verbose=False)["maxima"])
```

- **示例输出：**

```text
[2, 3, 8, 9]
```

### 14) `dnakit.epigenome.filter_summits_by_peakiness()`

- **作用：** 保留那些周围峰形足够宽、足够多样的峰顶，去掉过窄或过于平坦的峰。
- **处理方法：** 调用 MACS3 的 `SignalProcessing.enforce_peakyness`。峰的宽度至少 50 bp，且至少有 6 个不同取值。
- **API：** `dnakit.epigenome.filter_summits_by_peakiness(signal: 'Sequence[float]', summits: 'Sequence[int]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 信号列表与峰顶索引列表（从 0 开始，升序，不重复）。
- **输出：** `summits` 为保留的索引。
- **示例代码：**

```python
import math
from dnakit.epigenome import filter_summits_by_peakiness, find_signal_maxima

signal = [round(100 * math.exp(-((i - 60) ** 2) / 450), 4) for i in range(121)]
summits = find_signal_maxima(signal, window_size=21, verbose=False)["maxima"]
print(summits)
print(filter_summits_by_peakiness(signal, summits, verbose=False)["summits"])
```

- **示例输出：**

```text
[59]
[59]
```

### 15) `dnakit.epigenome.filter_summits_by_valleys()`

- **作用：** 检查相邻两个峰顶之间的谷值，保留谷足够低的峰顶。
- **处理方法：** 调用 MACS3 的 `SignalProcessing.enforce_valleys`。两峰之间的谷值必须不高于 `min_valley` 乘以较低峰顶的值。
- **API：** `dnakit.epigenome.filter_summits_by_valleys(signal: 'Sequence[float]', summits: 'Sequence[int]', *, min_valley: 'float' = 0.8, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 信号列表、峰顶索引列表，以及 `min_valley`（0 到 1）。
- **输出：** `summits` 为保留的索引。
- **示例代码：**

```python
from dnakit.epigenome import filter_summits_by_valleys

signal = [0, 5, 10, 5, 1, 4, 8, 4, 0]
print(filter_summits_by_valleys(signal, [2, 6], min_valley=0.8, verbose=False)["summits"])
print(filter_summits_by_valleys(signal, [2, 6], min_valley=0.1, verbose=False)["summits"])
```

- **示例输出：**

```text
[2, 6]
[2]
```

### 16) `dnakit.epigenome.refine_peak_summit()`

- **作用：** 在一个已确定的峰区间内，根据正负链 5' 端位置找出峰顶。
- **处理方法：** 调用 MACS3 的 `refinepeak_cmd.find_summit`，按窗口平滑后取读段末端最密集的位置。坐标与 MACS3 返回的一致。
- **API：** `dnakit.epigenome.refine_peak_summit(chrom: 'str', plus_positions: 'Sequence[int]', minus_positions: 'Sequence[int]', peak_start: 'int', peak_end: 'int', *, name: 'str' = 'peak', window_size: 'int' = 100, cutoff: 'float' = 5.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 染色体名、正负链 5' 端位置、峰的起止坐标、`window_size` 与 `cutoff`。
- **输出：** `summit_start`、`summit_end`、`name`、`score`。
- **示例代码：**

```python
from dnakit.epigenome import refine_peak_summit

result = refine_peak_summit(
    "chr1", [100, 110, 120], [105, 115], 90, 130, name="peak1", verbose=False
)
print(result["summit_start"], result["summit_end"], result["name"])
print(round(result["score"], 4))
```

- **示例输出：**

```text
101 102 peak1_F
0.8284
```

### 17) `dnakit.epigenome.max_duplicate_tags()`

- **作用：** 估计在给定 P 值下，同一位置上可能偶然出现的最多重复 tag 数，用于去重阈值。
- **处理方法：** 调用 MACS3 的 `callpeak_cmd.cal_max_dup_tags`，用二项分布的反向分布函数逐个试算。
- **API：** `dnakit.epigenome.max_duplicate_tags(genome_size: 'int', tag_count: 'int', *, p_value: 'float' = 1e-05, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 基因组大小、tag 总数、`p_value`。
- **输出：** `max_duplicate_tags` 为整数阈值。
- **示例代码：**

```python
from dnakit.epigenome import max_duplicate_tags

result = max_duplicate_tags(3_000_000_000, 10_000_000, verbose=False)
print(result["max_duplicate_tags"])
```

- **示例输出：**

```text
1
```

### 18) `dnakit.epigenome.guess_tag_file_format()`

- **作用：** 根据文件开头识别 MACS3 会用哪种读取器解析的 tag 文件格式。
- **处理方法：** 调用 MACS3 的 `IO.Parser.guess_parser`，返回第一个能识别文件的解析器。无法识别时报错。
- **API：** `dnakit.epigenome.guess_tag_file_format(path: 'str | Path', *, buffer_size: 'int' = 100000, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** tag 文件路径（BED、SAM 等）与读取缓冲大小 `buffer_size`。
- **输出：** `parser` 为解析器类名，`format` 为去掉 `Parser` 后缀的格式名。
- **示例代码：**

```python
from dnakit.epigenome import guess_tag_file_format

result = guess_tag_file_format("tags.bed", verbose=False)
print(result["parser"], result["format"])
```

- **示例输出：**

```text
BEDParser BED
```

### 19) `dnakit.epigenome.calculate_genotype_quality()`

- **作用：** 由三种基因型的自然对数似然计算 Phred 形式的基因型质量（GQ）。
- **处理方法：** 调用 MACS3 的 `VariantStat.calculate_GQ`，计算 `GQ = -10 log10((L2 + L3) / (L1 + L2 + L3))`，其中 `L = exp(lnL)`；MACS3 返回整数。
- **API：** `dnakit.epigenome.calculate_genotype_quality(log_likelihood_1: 'float', log_likelihood_2: 'float', log_likelihood_3: 'float', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 三个自然对数似然 `log_likelihood_1`、`log_likelihood_2`、`log_likelihood_3`。
- **输出：** `genotype_quality` 为整数。
- **示例代码：**

```python
from dnakit.epigenome import calculate_genotype_quality

result = calculate_genotype_quality(-10.0, -20.0, -30.0, verbose=False)
print(result["genotype_quality"])
```

- **示例输出：**

```text
43
```

### 20) `dnakit.epigenome.calculate_allele_model_likelihood()`

- **作用：** 用每条 read 的碱基质量计数，计算单个位点纯合或杂合模型的对数似然和 BIC。
- **处理方法：** 调用 MACS3 的 `VariantStat.CalModel_Homo`（纯合）或 `CalModel_Heter_noAS`（杂合，无等位偏好）。
- **API：** `dnakit.epigenome.calculate_allele_model_likelihood(top1_bq_t: 'Sequence[int]', top1_bq_c: 'Sequence[int]', top2_bq_t: 'Sequence[int]', top2_bq_c: 'Sequence[int]', *, model: 'str' = 'homozygous', verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 四个等长非负整数列表：第一与第二等位碱基的 T、C 质量；`model` 为 `"homozygous"` 或 `"heterozygous"`。
- **输出：** `log_likelihood`、`bic`、`reads`。
- **示例代码：**

```python
from dnakit.epigenome import calculate_allele_model_likelihood

result = calculate_allele_model_likelihood(
    [30, 25], [2, 3], [1, 1], [0, 2], model="heterozygous", verbose=False
)
print(result["reads"], round(result["log_likelihood"], 4), round(result["bic"], 4))
```

- **示例输出：**

```text
2 -6.5882 15.9489
```

### 21) `dnakit.epigenome.hmmratac_fragment_weights()`

- **作用：** 计算 HMMRATAC 中各片段长度属于 short、mono、di、tri 四类信号的权重。
- **处理方法：** 调用 MACS3 的 `HMMR_Signal_Processing.generate_weight_mapping`，四类各对应一个高斯分布（均值与标准差各 4 个值）。
- **API：** `dnakit.epigenome.hmmratac_fragment_weights(fragment_lengths: 'Sequence[int]', means: 'Sequence[float]', stddevs: 'Sequence[float]', *, min_fragment_probability: 'float' = 0.001, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 片段长度列表；四个均值 `means` 与四个标准差 `stddevs`，顺序为 short、mono、di、tri；`min_fragment_probability`。
- **输出：** `weights` 为每个片段长度一行，包含四类权重；`classes` 为类别名。
- **示例代码：**

```python
from dnakit.epigenome import hmmratac_fragment_weights

result = hmmratac_fragment_weights(
    [50, 150, 300],
    [50.0, 150.0, 300.0, 450.0],
    [20.0, 30.0, 40.0, 50.0],
    verbose=False,
)
for row in result["weights"]:
    print(row["fragment_length"], round(row["short"], 4), round(row["mono"], 4))
```

- **示例输出：**

```text
50 0.9974 0.0026
150 0.0 0.9993
300 0.0 0.0
```

### 22) `dnakit.epigenome.merge_pileup_arrays()`

- **作用：** 把两条覆盖度轨道按位置合并，每个位置取最大值、最小值或平均值。
- **处理方法：** 调用 MACS3 的 `Pileup.over_two_pv_array`，对两条 `[位置, 数值]` 轨道的并集逐点计算。
- **API：** `dnakit.epigenome.merge_pileup_arrays(first: 'Mapping[str, Sequence[float]]', second: 'Mapping[str, Sequence[float]]', *, reducer: 'str' = 'max', verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 两个映射，各含 `positions` 与 `values`；`reducer` 为 `"max"`、`"min"` 或 `"mean"`。
- **输出：** 合并后的 `positions` 与 `values`。
- **示例代码：**

```python
from dnakit.epigenome import merge_pileup_arrays

first = {"positions": [10, 20], "values": [1.0, 2.0]}
second = {"positions": [15, 20], "values": [3.0, 0.5]}
result = merge_pileup_arrays(first, second, reducer="max", verbose=False)
print(result["positions"])
print(result["values"])
```

- **示例输出：**

```text
[10, 15, 20]
[3.0, 3.0, 2.0]
```

## deepTools 函数（轨道、GC 与缩放因子）

### 23) `dnakit.epigenome.gc_fraction_from_twobit()`

- **作用：** 计算 2bit 基因组中一个区间的 GC 比例（或 G+C 碱基数）。
- **处理方法：** 用 py2bit 打开文件，调用 deepTools 的 `utilities.getGC_content`。N 碱基超过 5% 的区间会被拒绝。
- **API：** `dnakit.epigenome.gc_fraction_from_twobit(twobit_path: 'str | Path', chrom: 'str', start: 'int', end: 'int', *, fraction: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 2bit 文件、染色体名、0 为起点的半开区间 `[start, end)`；`fraction` 为真时返回比例。
- **输出：** `gc`（比例或碱基数）、`unit`、`start`、`end`。
- **示例代码：**

```python
import base64
from pathlib import Path
from dnakit.epigenome import gc_fraction_from_twobit

Path("genome.2bit").write_bytes(base64.b64decode("QydBGgAAAAABAAAAAAAAAARjaHIxGQAAABQAAAABAAAACgAAAAIAAAAAAAAAAAAAAJz1oPCc"))
print(gc_fraction_from_twobit("genome.2bit", "chr1", 0, 8, verbose=False)["gc"])
print(gc_fraction_from_twobit("genome.2bit", "chr1", 0, 8, fraction=False, verbose=False)["gc"])
```

- **示例输出：**

```text
0.75
6
```

### 24) `dnakit.epigenome.parse_deeptools_region()`

- **作用：** 解析 deepTools 风格的 `chrom:start:end[:tile]` 区域，并按染色体长度修正。
- **处理方法：** 调用 deepTools 的 `mapReduce.getUserRegion`。起点按瓦片大小对齐，终点截断到染色体长度。未知染色体会先被拒绝。
- **API：** `dnakit.epigenome.parse_deeptools_region(chrom_sizes: 'Mapping[str, int]', region: 'str', *, max_chunk_size: 'int' = 1000000, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `chrom_sizes`（染色体名到长度）、区域字符串、`max_chunk_size`。
- **输出：** `regions`、`start`、`end`、`chunk_size`。
- **示例代码：**

```python
from dnakit.epigenome import parse_deeptools_region

result = parse_deeptools_region({"chr2": 1000}, "chr2:10:1001", verbose=False)
print(result["regions"])
print(result["start"], result["end"], result["chunk_size"])
```

- **示例输出：**

```text
[{'chrom': 'chr2', 'length': 1000}]
10 1000 990
```

### 25) `dnakit.epigenome.estimate_bam_scale_factors()`

- **作用：** 由两份已建索引的 BAM 估计比例因子，用于让两个文库的覆盖度可比。
- **处理方法：** 调用 deepTools 的 `estimateScaleFactor`（SES 方法），在全基因组均匀取样 `number_of_samples` 个位点，比较读段计数。
- **API：** `dnakit.epigenome.estimate_bam_scale_factors(bam_paths: 'Sequence[str | Path]', *, bin_length: 'int' = 500, number_of_samples: 'int' = 1000, normalization_length: 'int' = 1000, avg_method: 'str' = 'median', verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 恰好两个 BAM 文件（各需 `.bai` 索引）；`bin_length`、`number_of_samples`、`normalization_length`、`avg_method`（`mean` 或 `median`）。
- **输出：** `size_factors`、`sites_sampled` 与其他估计值 `estimates`。
- **示例代码：**

```python
import random
import pysam
from pathlib import Path

random.seed(1)
sequence = "".join(random.choice("ACGT") for _ in range(400))
rows = []
for index in range(60):
    position = 1 + (index * 6) % 300
    bases = sequence[position - 1 : position + 19]
    rows.append(f"r{index}\t0\tchr1\t{position}\t60\t20M\t*\t0\t0\t{bases}\t{'I' * 20}")
sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:400\n" + "\n".join(rows) + "\n"
for name in ("scale_a", "scale_b"):
    Path(name + ".sam").write_text(sam)
    pysam.view("-b", "-o", name + ".unsorted.bam", name + ".sam", catch_stdout=False)
    pysam.sort("-o", name + ".bam", name + ".unsorted.bam")
    pysam.index(name + ".bam")
from dnakit.epigenome import estimate_bam_scale_factors

result = estimate_bam_scale_factors(
    ["scale_a.bam", "scale_b.bam"], bin_length=50, number_of_samples=40, verbose=False
)
print(result["size_factors"], result["sites_sampled"])
```

- **示例输出：**

```text
[1.0, 1.0] 36
```

### 26) `dnakit.epigenome.bedgraph_to_bigwig()`

- **作用：** 把已排序的 bedGraph 文件写成 bigWig 轨道文件。
- **处理方法：** 先检查染色体顺序、区间排序与不重叠、区间是否在染色体长度内，再调用 deepTools 的 `bedGraphToBigWig`（pyBigWig）。
- **API：** `dnakit.epigenome.bedgraph_to_bigwig(chrom_sizes: 'Mapping[str, int]', bedgraph_paths: 'Sequence[str | Path]', bigwig_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 染色体长度映射（顺序即 bigWig 头部顺序）、bedGraph 文件列表（按同样顺序）、新的 bigWig 输出路径。
- **输出：** `bigwig`、`chromosomes`、`intervals` 等。
- **示例代码：**

```python
from pathlib import Path
import pyBigWig
from dnakit.epigenome import bedgraph_to_bigwig

Path("signal.bw").unlink(missing_ok=True)
result = bedgraph_to_bigwig({"chr1": 100}, ["signal.bedgraph"], "signal.bw", verbose=False)
print(result["intervals"])
handle = pyBigWig.open("signal.bw")
print(list(handle.intervals("chr1")))
handle.close()
```

- **示例输出：**

```text
2
[(0, 10, 1.5), (10, 20, 2.0)]
```

### 27) `dnakit.epigenome.chop_regions_for_profile()`

- **作用：** 把外显子区间切分为不缩放的 5' 与 3' 两端以及缩放的主体，用于基因谱图（profile）。
- **处理方法：** 调用 deepTools 的 `heatmapper.chopRegions`。输入区间必须升序且不重叠。
- **API：** `dnakit.epigenome.chop_regions_for_profile(exons: 'Sequence[Sequence[int]]', *, left: 'int' = 0, right: 'int' = 0, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** `(start, end)` 区间列表；`left`、`right` 为两端各自的不缩放长度。
- **输出：** `left_bins`、`body_bins`、`right_bins`，以及 `pad_left`、`pad_right`（两端不足长度的补齐量）。
- **示例代码：**

```python
from dnakit.epigenome import chop_regions_for_profile

result = chop_regions_for_profile(
    [(0, 200), (300, 400), (800, 900)], left=50, right=50, verbose=False
)
print(result["left_bins"], result["right_bins"])
print(result["body_bins"])
```

- **示例输出：**

```text
[[0, 50]] [[850, 900]]
[[50, 200], [300, 400], [800, 850]]
```

### 28) `dnakit.epigenome.match_chromosome_name()`

- **作用：** 在名称列表中找到与请求名称对应的染色体名，允许增删 `chr` 前缀。
- **处理方法：** 调用 deepTools 的 `utilities.mungeChromosome`。找不到时返回 `None`。
- **API：** `dnakit.epigenome.match_chromosome_name(chrom: 'str', chromosome_names: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 请求的染色体名与可用名称列表。
- **输出：** `matched`：找到的名称或 `None`。
- **示例代码：**

```python
from dnakit.epigenome import match_chromosome_name

print(match_chromosome_name("chr1", ["1", "2"], verbose=False)["matched"])
print(match_chromosome_name("X", ["chr1"], verbose=False)["matched"])
```

- **示例输出：**

```text
1
None
```

### 29) `dnakit.epigenome.map_twobit_to_bam_chromosomes()`

- **作用：** 把 2bit 的染色体名映射到 BAM 的参考序列名，允许 `chr` 前缀不同。
- **处理方法：** 调用 deepTools 的 `utilities.tbitToBamChrName`。`complete` 表示所有 2bit 名称都得到了映射。
- **API：** `dnakit.epigenome.map_twobit_to_bam_chromosomes(twobit_names: 'Sequence[str]', bam_names: 'Sequence[str]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 2bit 染色体名列表与 BAM 参考序列名列表。
- **输出：** `mapping`（2bit 名到 BAM 名）与 `complete`。
- **示例代码：**

```python
from dnakit.epigenome import map_twobit_to_bam_chromosomes

result = map_twobit_to_bam_chromosomes(["chr1", "chr2"], ["1", "2"], verbose=False)
print(result["mapping"], result["complete"])
partial = map_twobit_to_bam_chromosomes(["chr1", "chrX"], ["chr1", "chr2"], verbose=False)
print(partial["complete"])
```

- **示例输出：**

```text
{'chr1': '1', 'chr2': '2'} True
False
```

### 30) `dnakit.epigenome.fingerprint_jensen_shannon()`

- **作用：** 计算 ChIP 与对照（input）覆盖度的 Jensen-Shannon 指标，即 plotFingerprint 报告的 JSD。
- **处理方法：** 调用 deepTools 的 `plotFingerprint.getJSDcommon`，比较两者的累积覆盖曲线。计数按输入顺序使用。
- **API：** `dnakit.epigenome.fingerprint_jensen_shannon(chip_counts: 'Sequence[float]', input_counts: 'Sequence[float]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** ChIP 与 input 的每个 bin 计数列表（等长，至少 2 个）。
- **输出：** `jsd` 为数值，`bins` 为 bin 数。
- **示例代码：**

```python
from dnakit.epigenome import fingerprint_jensen_shannon

result = fingerprint_jensen_shannon([0, 1, 2, 3, 10], [1, 1, 1, 1, 1], verbose=False)
print(round(result["jsd"], 6), result["bins"])
print(fingerprint_jensen_shannon([1, 2, 3, 4], [1, 2, 3, 4], verbose=False)["jsd"])
```

- **示例输出：**

```text
0.283567 5
0.0
```

### 31) `dnakit.epigenome.estimate_count_size_factors()`

- **作用：** 由 bin×样本计数矩阵估计大小因子，结果可直接用作 bamCoverage 的比例因子。
- **处理方法：** 调用 deepTools 的 `countReadsPerBin.estimateSizeFactors`，即 DESeq2 方法（只使用无零计数的行），返回其倒数。
- **API：** `dnakit.epigenome.estimate_count_size_factors(counts: 'Sequence[Sequence[float]]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 等宽的非负数值矩阵（每行一个 bin，每列一个样本，至少 2 列）。
- **输出：** `size_factors`，每个样本一个值。
- **示例代码：**

```python
from dnakit.epigenome import estimate_count_size_factors

matrix = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [0, 10, 0], [10, 5, 100]]
result = estimate_count_size_factors(matrix, verbose=False)
print([round(value, 4) for value in result["size_factors"]])
```

- **示例输出：**

```text
[1.305, 0.9932, 0.783]
```

### 32) `dnakit.epigenome.bigwig_chromosome_sizes()`

- **作用：** 读取 bigWig 文件中的染色体名与长度，并列出各文件之间不一致的条目。
- **处理方法：** 调用 deepTools 的 `getScorePerBigWigBin.getChromSizes`，通过 pyBigWig 读取每个文件的染色体表。
- **API：** `dnakit.epigenome.bigwig_chromosome_sizes(bigwig_paths: 'Sequence[str | Path]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **输入：** 一到 32 个 bigWig 文件路径。
- **输出：** `chromosomes`（共同的染色体）与 `non_matching`（不一致的条目）。
- **示例代码：**

```python
import pyBigWig
from pathlib import Path
from dnakit.epigenome import bigwig_chromosome_sizes

Path("track.bw").unlink(missing_ok=True)
handle = pyBigWig.open("track.bw", "w")
handle.addHeader([("chr1", 100)])
handle.addEntries(["chr1"], [0], ends=[10], values=[1.0])
handle.close()
print(bigwig_chromosome_sizes(["track.bw"], verbose=False)["chromosomes"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'length': 100}]
```
