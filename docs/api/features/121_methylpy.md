# methylpy：亚硫酸氢盐比对、甲基化位点检出与差异甲基化区域查找

这一组函数是 `dnakit.epigenome` 中基于 [methylpy](https://github.com/yupenghe/methylpy) 的适配器，methylpy 是用于全基因组亚硫酸氢盐测序分析的 Python 流程。它们构建 C 转 T 与 G 转 A 的转换参考序列，比对单端和双端亚硫酸氢盐测序读段，合并并去重比对结果，修剪接头，把甲基化胞嘧啶检出到 allc 表中，过滤并压缩这些表，把它们转换为 bigWig 轨道，并查找差异甲基化区域（DMR）。每个函数调用对应的 methylpy 函数，只转发你设置的参数，其余沿用 methylpy 自己的默认值。DNAKit 额外做了三件事：先检查外部程序（缺少程序时抛出带安装提示的 `BackendUnavailableError`，而不是让 methylpy 调用 `sys.exit`），把 methylpy 的失败抛为 `BackendExecutionError`，并返回可 JSON 序列化的字典。其中许多函数封装了运行 `samtools`、`bowtie2`、`bgzip` 等程序的比对与索引步骤，并会在输出前缀旁写出多个文件。

methylpy 1.4.7 的包元数据附带 Apache License 2.0，因此该 Python 包没有 copyleft 方面的限制。DMR 函数使用的内置 `run_rms_tests.out` 二进制文件链接 GNU Scientific Library（`libgsl.so.0`）。

依赖都是可选的，只有调用对应函数时才会导入。用 `pip install "dnakit[bio-methylpy]"` 安装 methylpy，用 `pip install "dnakit[bio-cutadapt]"` 安装 `cutadapt`；`samtools`、`bowtie2` 或 `bowtie`、`bgzip`、`tabix`、`wigToBigWig`、Picard 与 GSL 是需要另行安装的外部程序和库。示例从 `dnakit.epigenome` 导入。

本环境装有 methylpy 与 cutadapt，但没有 `samtools`、`bowtie2`、`bgzip`、`tabix`、`wigToBigWig` 和 GSL。因此有五个示例执行真实计算（`methylpy_filter_allc_files`、两个 `find_multi_mappers` 函数和两个 `quality_trim` 函数）。其余示例展示缺少程序或输入时函数的真实行为：它们把程序路径指向一个不存在的文件夹（`path_to_samtools="no_such_dir"`），因此在任何机器上输出都一样；每个这样的示例之后都有一个不执行的“安装依赖后的用法”代码块，给出真实可用的调用。

## 示例数据

示例使用三个很小的声明文件：一个五行的 allc 表、一个单序列参考 FASTA 和两条 FASTQ 读段（第一条带 3' 接头）。其他输入由示例自己创建，输出写入临时目录，或根本不会生成。

**示例文件：** `allc.tsv`

```text
chr1	100	+	CGG	10	12	1
chr1	105	-	CHH	2	20	0
chr1	300	+	CAG	3	4	1
chr2	50	+	CGT	0	8	0
chr2	70	-	CCC	5	6	1
```

**示例文件：** `ref.fa`

```text
>chr1
ACGTACGTTTGGAACCGGTTACGTACGTTTGGAACCGGTT
```

**示例文件：** `reads.fq`

```text
@r1
ACGTACGTTTGGAAGATCGGAAGAGCACACGTCTGAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@r2
TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## 1) `dnakit.epigenome.methylpy_bgzip_allc_file()`

- **作用：** 用 bgzip 压缩 allc 甲基化表，并用 tabix 建立索引。
- **处理方法：** 检查文件存在并检查 `bgzip` 与 `tabix`，然后调用 `methylpy.utilities.bgzip_allc_file`。上游写出分块压缩的副本 `<allc_file>.bgz`（输入已以 `.gz` 结尾时为用 `.bgz` 替换 `.gz` 的名称），并用 `tabix -s 1 -b 2 -e 2` 建立索引，生成 `.tbi` 文件。输入文件会保留。
- **API：** `dnakit.epigenome.methylpy_bgzip_allc_file(allc_file: 'str', *, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, buffer_line_number: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `bgzip` 与 `tabix`（htslib）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `allc_file`：已存在的 allc 表。`path_to_bgzip` 与 `path_to_tabix` 是程序所在的文件夹（methylpy 会追加程序名）。`buffer_line_number` 设置每次读取的行数。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_bgzip_allc_file
from dnakit.exceptions import BackendUnavailableError

try:
    methylpy_bgzip_allc_file("allc.tsv", path_to_bgzip="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - bgzip was not found at no_such_dir/bgzip.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_bgzip_allc_file

result = methylpy_bgzip_allc_file("allc_sample_chr1.tsv", verbose=False)
print(result["status"])  # allc_sample_chr1.tsv.bgz and allc_sample_chr1.tsv.bgz.tbi are written next to the input
```

## 2) `dnakit.epigenome.methylpy_convert_allc_to_bigwig()`

- **作用：** 把 allc 表转换为按窗口统计的甲基化 bigWig 轨道。
- **处理方法：** 检查输入文件，检查 `samtools` 与 `wigToBigWig`，然后调用 `methylpy.utilities.convert_allc_to_bigwig`，当 `<reference_fasta>.fai` 不存在时，上游用 `samtools faidx` 为 FASTA 建索引，写出临时的 `<output_file>.chrom_size` 与 `<output_file>.wig`，再运行 `wigToBigWig`。在 wig 文件中，每个 `bin_size` 碱基（默认 100）的窗口保存所选上下文的甲基化水平，即甲基化读段数之和除以覆盖读段数之和，临时文件随后被删除。
- **API：** `dnakit.epigenome.methylpy_convert_allc_to_bigwig(input_allc_file: 'str', output_file: 'str', reference_fasta: 'str', *, mc_type: 'str | None' = None, bin_size: 'int | None' = None, path_to_wigToBigWig: 'str | None' = None, path_to_samtools: 'str | None' = None, min_bin_sites: 'int | None' = None, min_bin_cov: 'int | None' = None, max_site_cov: 'int | None' = None, min_site_cov: 'int | None' = None, remove_chr_prefix: 'bool | None' = None, add_chr_prefix: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools` 与 `wigToBigWig`（UCSC）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `input_allc_file` 与 `reference_fasta` 必须存在；`output_file` 是 bigWig 路径。可选：`mc_type`、`bin_size`、`min_bin_sites`、`min_bin_cov`、`min_site_cov`、`max_site_cov`、`remove_chr_prefix`、`add_chr_prefix`、`path_to_wigToBigWig`、`path_to_samtools`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（调用后存在的已声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_convert_allc_to_bigwig
from dnakit.exceptions import BackendUnavailableError

try:
    methylpy_convert_allc_to_bigwig("allc.tsv", "cg.bw", "ref.fa", path_to_samtools="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_convert_allc_to_bigwig

result = methylpy_convert_allc_to_bigwig(
    "allc_sample_chr1.tsv", "sample_cg.bw", "genome.fa", mc_type="CGN", bin_size=100, verbose=False
)
print(result["outputs"])
```

## 3) `dnakit.epigenome.methylpy_dmrfind()`

- **作用：** 根据多个 allc 文件查找样本之间的差异甲基化区域。
- **处理方法：** 检查内置的 `run_rms_tests.out` 二进制文件存在，检查每个 allc 文件存在，然后调用 `methylpy.DMRfind.DMRfind`。上游用均方根检验寻找差异甲基化位点，合并 `dmr_max_dist` 以内的相邻位点，并通过内置二进制文件运行置换检验。只给一个文件时上游会以“Needs more than one allc file as input”退出，DNAKit 将其报告为 `BackendExecutionError`。
- **API：** `dnakit.epigenome.methylpy_dmrfind(allc_files: 'Sequence[str]', output_prefix: 'str', *, samples: 'Sequence[str] | None' = None, mc_type: 'str | None' = None, chroms: 'Sequence[str] | None' = None, num_procs: 'int | None' = None, min_cov: 'int | None' = None, dmr_max_dist: 'int | None' = None, sig_cutoff: 'float | None' = None, num_sims: 'int | None' = None, num_sig_tests: 'int | None' = None, seed: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`）；内置的 `run_rms_tests.out` 二进制文件需要 GSL 共享库 `libgsl.so.0` 和 `libgslcblas.so.0`（`apt install libgsl-dev` 或 `conda install -c conda-forge gsl`）。
- **输入：** `allc_files`：至少两个已存在的 allc 表的列表；`output_prefix`：结果文件前缀。可选：`samples`（样本名，每个文件一个）、`mc_type`、`chroms`、`num_procs`、`min_cov`、`dmr_max_dist`、`sig_cutoff`、`num_sims`、`num_sig_tests`、`seed`、`keep_temp_files`。`mc_type` 原样传给 methylpy（其默认值为 `CGN`）。
- **输出：** 上述结果字典（`outputs` 为空列表；methylpy 会写出 `<output_prefix>_DMRs.tsv` 等文件）。输入文件不存在时抛出 `ConfigurationError`。
- **示例代码：**

```python
from dnakit.epigenome import methylpy_dmrfind
from dnakit.exceptions import ConfigurationError

try:
    methylpy_dmrfind(["missing_a.tsv", "missing_b.tsv"], "dmr", samples=["a", "b"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "-", error.message.split(":")[0])
```

- **示例输出：**

```text
ConfigurationError - allc_files must be an existing file
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_dmrfind

result = methylpy_dmrfind(
    ["allc_a.tsv.gz", "allc_b.tsv.gz"], "dmr", samples=["a", "b"], mc_type="CGN", num_procs=2, verbose=False
)
print(result["status"])  # dmr_DMRs.tsv is written
```

## 4) `dnakit.epigenome.methylpy_pe_call_methylated_sites()`

- **作用：** 从排序后的双端 BAM 文件检出甲基化胞嘧啶。
- **处理方法：** 检查 BAM 与参考序列存在，检查 `samtools`（`bgzip=True` 时还检查 `bgzip` 与 `tabix`），然后调用 `methylpy.call_mc_pe.call_methylated_sites_pe`。上游按染色体生成 pileup，用 `unmethylated_control` 估计未转化率，并以 `sig_cutoff` 调整后的 p 值做二项检验，每条染色体写出一个 `allc_<sample>_<chrom>.tsv`。
- **API：** `dnakit.epigenome.methylpy_pe_call_methylated_sites(inputf: 'str', sample: 'str', reference_fasta: 'str', *, unmethylated_control: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_mpileup_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, buffer_line_number: 'int | None' = None, min_mapq: 'int | None' = None, min_cov: 'int | None' = None, binom_test: 'bool | None' = None, path_to_samtools: 'str | None' = None, remove_chr_prefix: 'bool | None' = None, add_snp_info: 'bool | None' = None, sort_mem: 'str | None' = None, path_to_files: 'str | None' = None, min_base_quality: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（`bgzip=True` 时还有 `bgzip` 与 `tabix`）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `inputf`：排序后的 BAM；`sample`：输出文件中使用的名称；`reference_fasta`：已存在的 FASTA。可选：`unmethylated_control`（`chrom:start-end` 区域或一个数值）、`sig_cutoff`、`num_procs`、`num_upstr_bases`、`num_downstr_bases`、`min_mapq`、`min_cov`、`min_base_quality`、`binom_test`、`generate_mpileup_file`、`compress_output`、`bgzip`、`remove_chr_prefix`、`add_snp_info`、`sort_mem`、`buffer_line_number`、`path_to_files`、`path_to_samtools`、`path_to_bgzip`、`path_to_tabix`、`keep_temp_files`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_call_methylated_sites
from dnakit.exceptions import BackendUnavailableError

Path("reads_pe.bam").touch()
try:
    methylpy_pe_call_methylated_sites(
        "reads_pe.bam", "s1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_call_methylated_sites

result = methylpy_pe_call_methylated_sites(
    "s1_processed_reads_no_clonal.bam", "s1", "genome.fa",
    unmethylated_control="chrC:", num_procs=4, verbose=False,
)
print(result["status"])  # allc_s1_<chrom>.tsv.gz files are written
```

## 5) `dnakit.epigenome.methylpy_pe_find_multi_mappers()`

- **作用：** 从双端比对中筛出正确配对且质量足够的读段，并还原原始读段序列。
- **处理方法：** 调用 `methylpy.call_mc_pe.find_multi_mappers_pe`。它读取比对程序写出的 SAM 文件，跳过头部行、未正确配对的读段（标志位 2）、比对质量低于 `max(3, min_mapq)` 的读段以及补充比对，并利用读段名中 `!` 之后编码的位置还原亚硫酸氢盐转换前的序列。保留的行被轮流写入 `<output>_sorted_<i>`（`i < num_procs`）。除非 `keep_temp_files=True`，输入 SAM 文件会被删除。
- **API：** `dnakit.epigenome.methylpy_pe_find_multi_mappers(inputf: 'str', output: 'str', *, num_procs: 'int | None' = None, min_mapq: 'int | None' = None, keep_temp_files: 'bool | None' = None, append: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 只需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），不运行外部程序。输入必须是读段名以 `!<编码位置>` 结尾的 SAM 文件，即 methylpy 比对步骤写出的格式。
- **输入：** `inputf`：已存在的 SAM 文件；`output`：结果文件前缀。可选：`num_procs`（输出文件数）、`min_mapq`、`keep_temp_files`、`append`（追加到已有输出文件而不是覆盖）。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（此处恒为空列表，因为该函数没有声明单个输出路径：methylpy 按前缀命名并写出文件）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_pe_find_multi_mappers

folder = Path(tempfile.mkdtemp())
sam = folder / "pairs.sam"
sam.write_text(
    "@HD\tVN:1.0\n"
    "p1!\t99\tchr1\t10\t42\t8M\t=\t30\t28\tACGTACGT\tIIIIIIII\n"
    "p1!\t147\tchr1\t30\t42\t8M\t=\t10\t-28\tTTGGTTAA\tIIIIIIII\n"
    "p2!\t65\tchr1\t10\t42\t8M\t=\t30\t28\tACGTACGT\tIIIIIIII\n"
)
result = methylpy_pe_find_multi_mappers(
    str(sam), str(folder / "pe"), num_procs=2, keep_temp_files=True, verbose=False
)
print(result["status"], result["outputs"])
for index in range(2):
    rows = (folder / f"pe_sorted_{index}").read_text().splitlines()
    print(index, [row.split("\t")[:2] for row in rows])
```

- **示例输出：**

```text
completed []
0 [['p1', '99']]
1 [['p1', '147']]
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_find_multi_mappers

result = methylpy_pe_find_multi_mappers("s1_forward.sam", "s1_forward", num_procs=4, verbose=False)
# s1_forward_sorted_0 ... s1_forward_sorted_3 are written; s1_forward.sam is deleted
```

## 6) `dnakit.epigenome.methylpy_pe_flip_read2_strand()`

- **作用：** 翻转 BAM 文件中每个已比对读段对里 read 2 的链方向。
- **处理方法：** 检查输入存在且 `samtools` 可用，然后调用 `methylpy.call_mc_pe.flip_read2_strand`，重写 BAM，使一对读段的两条都报告原始链。
- **API：** `dnakit.epigenome.methylpy_pe_flip_read2_strand(input_file: 'str', output_file: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `input_file`：已存在的双端 BAM；`output_file`：要写出的 BAM。`path_to_samtools` 是 samtools 所在的文件夹。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（调用后存在的已声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_flip_read2_strand
from dnakit.exceptions import BackendUnavailableError

Path("pairs.bam").touch()
try:
    methylpy_pe_flip_read2_strand("pairs.bam", "flipped.bam", path_to_samtools="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_flip_read2_strand

result = methylpy_pe_flip_read2_strand("s1_pe.bam", "s1_pe_flipped.bam", verbose=False)
print(result["outputs"])
```

## 7) `dnakit.epigenome.methylpy_pe_merge_bam_files()`

- **作用：** 把多个双端 BAM 文件合并为一个，并修正头部。
- **处理方法：** 检查输入存在且 `samtools` 可用，然后调用 `methylpy.call_mc_pe.merge_bam_files`，它在当前目录写出临时的 `header.sam`（第一个文件的头部加每个输入一行 `@RG`，随后删除），并运行 `samtools merge -r -h header.sam`。每个输入文件名都必须包含 `.bam`。
- **API：** `dnakit.epigenome.methylpy_pe_merge_bam_files(input_files: 'Sequence[str]', output: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `input_files`：非空的已存在 BAM 文件列表；`output`：合并后的 BAM；`path_to_samtools`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（调用后存在的已声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_merge_bam_files
from dnakit.exceptions import BackendUnavailableError

Path("lib1.bam").touch()
Path("lib2.bam").touch()
try:
    methylpy_pe_merge_bam_files(
        ["lib1.bam", "lib2.bam"], "merged.bam", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_merge_bam_files

result = methylpy_pe_merge_bam_files(
    ["s1_lib1_no_clonal.bam", "s1_lib2_no_clonal.bam"], "s1_merged.bam", verbose=False
)
print(result["outputs"])
```

## 8) `dnakit.epigenome.methylpy_pe_merge_sorted_multimap()`

- **作用：** 合并一个文库的排序后双端比对文件，只保留唯一比对的读段对。
- **处理方法：** 检查文件与参考序列存在且 `samtools` 可用，然后调用 `methylpy.call_mc_pe.merge_sorted_multimap_pe`。当 `<reference_fasta>.fai` 不存在时，上游用 `samtools faidx` 为参考序列建索引，根据索引写出 SAM 头部，逐条合并排序后的文件并只保留在所有文件中仅出现一次的读段对，最后转换为 `<prefix>_processed_reads.bam`。
- **API：** `dnakit.epigenome.methylpy_pe_merge_sorted_multimap(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `current_library`：文库 ID；`files`：`methylpy_pe_find_multi_mappers` 写出的 `_sorted_<i>` 文件；`prefix`：输出前缀；`reference_fasta`：已存在的 FASTA。`path_to_samtools` 可选。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_merge_sorted_multimap
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_pe_merge_sorted_multimap(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_merge_sorted_multimap

result = methylpy_pe_merge_sorted_multimap(
    "lib1", ["s1_lib1_sorted_0", "s1_lib1_sorted_1"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 9) `dnakit.epigenome.methylpy_pe_merge_sorted_multimap_max_mapq()`

- **作用：** 合并一个文库的排序后双端比对文件，并为多重比对读段保留 MAPQ 最高的读段对。
- **处理方法：** 检查文件与参考序列存在且 `samtools` 可用，然后调用 `methylpy.call_mc_pe.merge_sorted_multimap_pe_max_mapq`。当 `<reference_fasta>.fai` 不存在时，上游用 `samtools faidx` 为参考序列建索引，根据索引写出 SAM 头部，逐条合并排序后的文件并对出现多次的读段保留比对质量最高的读段对，最后转换为 `<prefix>_processed_reads.bam`。
- **API：** `dnakit.epigenome.methylpy_pe_merge_sorted_multimap_max_mapq(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** 与 `methylpy_pe_merge_sorted_multimap` 相同：`current_library`、`files`、`prefix`、已存在的 `reference_fasta`，以及可选的 `path_to_samtools`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_merge_sorted_multimap_max_mapq
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_pe_merge_sorted_multimap_max_mapq(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_merge_sorted_multimap_max_mapq

result = methylpy_pe_merge_sorted_multimap_max_mapq(
    "lib1", ["s1_lib1_sorted_0"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 10) `dnakit.epigenome.methylpy_pe_run_alignment()`

- **作用：** 把一个文库的双端读段比对到两个转换参考序列上，并合并两条链的结果。
- **处理方法：** 检查文件与程序，然后调用 `methylpy.call_mc_pe.run_alignment_pe`。上游用所选比对程序（默认 `bowtie2`，也可用 `bowtie`）把读段分别比对到 C 转 T 参考序列和 G 转 A 参考序列，合并两条链并保留只在一处比对上的读段。
- **API：** `dnakit.epigenome.methylpy_pe_run_alignment(current_library: 'str', library_read1_files: 'Sequence[str]', library_read2_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, path_to_aligner: 'str | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, num_procs: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools` 与 `bowtie2`（`aligner="bowtie"` 时为 `bowtie`）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `current_library`：文库 ID；`library_read1_files` 与 `library_read2_files`：已存在的 FASTQ 文件列表；`sample`；`forward_reference` 与 `reverse_reference`：`methylpy_se_build_ref` 构建的转换 FASTA（必须存在）；`reference_fasta`。可选：`path_to_output`、`aligner`、`aligner_options`、`path_to_aligner`、`path_to_samtools`、`merge_by_max_mapq`、`min_mapq`、`num_procs`、`keep_temp_files`、`sort_mem`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_run_alignment
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_pe_run_alignment(
        "lib1", ["reads.fq"], ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_run_alignment

result = methylpy_pe_run_alignment(
    "lib1", ["lib1_R1.fastq.gz"], ["lib1_R2.fastq.gz"], "s1",
    "genome_f.fasta", "genome_r.fasta", "genome.fa", num_procs=4, verbose=False,
)
```

## 11) `dnakit.epigenome.methylpy_pe_run_mapping()`

- **作用：** 对一个文库的双端亚硫酸氢盐读段完成修剪、比对与合并。
- **处理方法：** 检查文件与程序，然后调用 `methylpy.call_mc_pe.run_mapping_pe`。上游用所选比对程序（默认 `bowtie2`，也可用 `bowtie`）把读段分别比对到 C 转 T 参考序列和 G 转 A 参考序列，在可选的 cutadapt 接头修剪（`trim_reads`，默认开启）之后，合并两条链并保留只在一处比对上的读段。
- **API：** `dnakit.epigenome.methylpy_pe_run_mapping(current_library: 'str', library_read1_files: 'Sequence[str]', library_read2_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, num_procs: 'int | None' = None, trim_reads: 'bool | None' = None, path_to_cutadapt: 'str | None' = None, adapter_seq_read1: 'str | None' = None, adapter_seq_read2: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`、`bowtie2`（或 `bowtie`）与 `cutadapt`（`trim_reads=False` 时不需要）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** 与 `methylpy_pe_run_alignment` 相同（参考文件必须存在），另有修剪选项 `trim_reads`、`path_to_cutadapt`、`adapter_seq_read1`、`adapter_seq_read2`、`min_qual_score` 与 `min_read_len`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_run_mapping
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_pe_run_mapping(
        "lib1", ["reads.fq"], ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_run_mapping

result = methylpy_pe_run_mapping(
    "lib1", ["lib1_R1.fastq.gz"], ["lib1_R2.fastq.gz"], "s1",
    "genome_f.fasta", "genome_r.fasta", "genome.fa", num_procs=4, verbose=False,
)
```

## 12) `dnakit.epigenome.methylpy_pe_run_methylation_pipeline()`

- **作用：** 运行从 FASTQ 文件到 allc 表的完整双端亚硫酸氢盐流程。
- **处理方法：** 检查输入与程序，然后调用 `methylpy.call_mc_pe.run_methylation_pipeline_pe`。上游依次修剪接头，把每个文库比对到两个转换参考序列，合并文库，可选地用 Picard 去除克隆读段，并按染色体检出甲基化位点。这是串联本页其他双端函数的入口。
- **API：** `dnakit.epigenome.methylpy_pe_run_methylation_pipeline(read1_files: 'Sequence[str]', read2_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, libraries: 'Sequence[str] | None' = None, unmethylated_control: 'str | None' = None, path_to_output: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, sort_mem: 'str | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_allc_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, binom_test: 'bool | None' = None, min_cov: 'int | None' = None, trim_reads: 'bool | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, remove_clonal: 'bool | None' = None, remove_chr_prefix: 'bool | None' = None, pbat: 'bool | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, path_to_cutadapt: 'str | None' = None, path_to_picard: 'str | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, adapter_seq_read1: 'str | None' = None, adapter_seq_read2: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`、`bowtie2`（或 `bowtie`）、`cutadapt`（`trim_reads=False` 时不需要）、`bgzip` 与 `tabix`（仅 `bgzip=True` 时）；`remove_clonal=True` 还需要 Picard 与 Java，它们不会被预先检查（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `read1_files` 与 `read2_files`：顺序一一对应的已存在 FASTQ 文件列表；`sample`；`forward_reference` 与 `reverse_reference`（已存在的转换 FASTA）；`reference_fasta`。可选：`libraries`、`unmethylated_control`、`path_to_output`、`sig_cutoff`、`num_procs`、`min_cov`、`trim_reads`、`aligner`、`aligner_options`、`merge_by_max_mapq`、`min_mapq`、`remove_clonal`、`remove_chr_prefix`、`pbat`、`generate_allc_file`、`compress_output`、`bgzip`、`binom_test`、各 `path_to_...` 参数、接头与质量选项以及 `keep_temp_files`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_pe_run_methylation_pipeline
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_pe_run_methylation_pipeline(
        ["reads.fq"], ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_pe_run_methylation_pipeline

result = methylpy_pe_run_methylation_pipeline(
    ["lib1_R1.fastq.gz"], ["lib1_R2.fastq.gz"], "s1",
    "genome_f.fasta", "genome_r.fasta", "genome.fa",
    unmethylated_control="chrC:", num_procs=8, verbose=False,
)
```

## 13) `dnakit.epigenome.methylpy_quality_trim_pe()`

- **作用：** 用 cutadapt 修剪双端亚硫酸氢盐读段的接头和低质量末端。
- **处理方法：** 检查输入存在且输出尚不存在，检查 `cutadapt`，然后用单元素列表调用 `methylpy.call_mc_pe.quality_trim_pe`。上游会错误处理普通字符串路径（把列表转换赋给了错误的变量名），所以 DNAKit 总是传入列表。上游用文件名拼出 shell 命令行，因此路径中应避免特殊字符。默认接头是 Illumina 序列 `AGATCGGAAGAGCACACGTCTGAAC`（read 1）与 `AGATCGGAAGAGCGTCGTGTAGGGA`（read 2）；`min_qual_score` 与 `min_read_len` 的默认值为 10 与 30。示例中两个输出来自同一个输入文件，所以只有 read 1 被去掉接头（文件中带的是 read 1 的接头，而不是 read 2 的接头）。
- **API：** `dnakit.epigenome.methylpy_quality_trim_pe(inputf_read1: 'str', outputf_read1: 'str', inputf_read2: 'str', outputf_read2: 'str', *, quality_base: 'int | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, adapter_seq_read1: 'str | None' = None, adapter_seq_read2: 'str | None' = None, num_procs: 'int | None' = None, input_format: 'str | None' = None, error_rate: 'float | None' = None, max_adapter_removal: 'int | None' = None, overlap_length: 'int | None' = None, zero_cap: 'bool | None' = None, path_to_cutadapt: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`）和外部程序 `cutadapt`（额外依赖 `dnakit[bio-cutadapt]`）。上游通过 `path_to_cutadapt`（文件夹）或 PATH 找到它；示例传入当前 Python 所在的文件夹，pip 会把它装在那里。
- **输入：** `inputf_read1`、`inputf_read2`：已存在的 FASTQ 文件；`outputf_read1`、`outputf_read2`：必须尚不存在的新路径。可选：`quality_base`、`min_qual_score`、`min_read_len`、`adapter_seq_read1`、`adapter_seq_read2`、`num_procs`、`input_format`、`error_rate`、`max_adapter_removal`、`overlap_length`、`zero_cap`、`path_to_cutadapt`。
- **输出：** 上述结果字典，其中 `outputs` 是两个修剪后的 FASTQ 路径。
- **示例代码：**

```python
import sys
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_quality_trim_pe

folder = Path(tempfile.mkdtemp())
result = methylpy_quality_trim_pe(
    "reads.fq", str(folder / "r1.fq"), "reads.fq", str(folder / "r2.fq"),
    min_read_len=10,
    path_to_cutadapt=str(Path(sys.executable).parent),
    verbose=False,
)
print(result["status"], len(result["outputs"]))
for name in ("r1.fq", "r2.fq"):
    lines = (folder / name).read_text().splitlines()
    print(name, [lines[i + 1] for i in range(0, len(lines), 4)])
```

- **示例输出：**

```text
completed 2
r1.fq ['ACGTACGTTTGGA', 'TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT']
r2.fq ['ACGTACGTTTGGAAGATCGGAAGAGCACACGTCTGAAAA', 'TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT']
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_quality_trim_pe

result = methylpy_quality_trim_pe(
    "lib1_R1.fastq", "lib1_R1_trimmed.fastq", "lib1_R2.fastq", "lib1_R2_trimmed.fastq",
    min_qual_score=20, min_read_len=30, verbose=False,
)
```

## 14) `dnakit.epigenome.methylpy_quality_trim_se()`

- **作用：** 用 cutadapt 修剪单端亚硫酸氢盐读段的接头和低质量末端。
- **处理方法：** 检查输入存在、输出尚不存在、格式和接头，找到 `cutadapt`，并对每个文件以参数列表（不经 shell）运行一次（这取代了 methylpy 用插值字符串拼 shell 命令的 `quality_trim`）。命令为 `cutadapt --quiet [-z] [-f F] [-e E] [-n N] [-O O] [-m MIN_LEN] [-q MIN_QUAL] [--quality-base=B] -a ADAPTER -o OUT IN`。只有取值为真时才添加对应选项。与另一个修剪函数不同，`min_qual_score` 与 `min_read_len` 不使用 methylpy 的默认值：除非你设置，否则沿用 cutadapt 自己的默认行为（不修剪质量，不过滤长度）。
- **API：** `dnakit.epigenome.methylpy_quality_trim_se(inputf: 'Sequence[str]', output: 'Sequence[str]', *, quality_base: 'int | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, adapter_seq: 'str' = 'AGATCGGAAGAGCACACGTCTG', input_format: 'str | None' = None, error_rate: 'float | None' = None, max_adapter_removal: 'int | None' = None, overlap_length: 'int | None' = None, zero_cap: 'bool' = False, path_to_cutadapt: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `cutadapt`（额外依赖 `dnakit[bio-cutadapt]`），它应在 PATH 中、在当前 Python 所在的文件夹中，或通过文件夹 `path_to_cutadapt` 指定。它本身不需要 methylpy，只是用来报告 methylpy 的版本。
- **输入：** `inputf`：已存在的 FASTQ 或 FASTA 文件列表；`output`：等长的新路径列表，这些路径必须尚不存在。`adapter_seq`（默认 `AGATCGGAAGAGCACACGTCTG`，只允许 IUPAC 字母）、`input_format`（`fasta`、`fastq` 或 `sra-fastq`）、`min_read_len`、`min_qual_score`、`quality_base`、`error_rate`、`max_adapter_removal`、`overlap_length`、`zero_cap`、`path_to_cutadapt`。
- **输出：** 字典，包含 `status`、`backend`、`version`、`function`、`outputs`（存在的修剪后文件）和 `log_tail`（每个文件 cutadapt stderr 的最后 4000 个字符）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_quality_trim_se

output = Path(tempfile.mkdtemp()) / "trimmed.fq"
result = methylpy_quality_trim_se(["reads.fq"], [output], min_read_len=10, verbose=False)
lines = output.read_text().splitlines()
print(result["status"], len(result["outputs"]), len(lines) // 4)
for index in range(0, len(lines), 4):
    print(lines[index], lines[index + 1])
```

- **示例输出：**

```text
completed 1 2
@r1 ACGTACGTTTGGA
@r2 TTGGTTGGTTAAGGCCAAGGTTTTAAGGCCTTAAGGAATT
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_quality_trim_se

result = methylpy_quality_trim_se(
    ["lib1.fastq", "lib2.fastq"], ["lib1_trimmed.fastq", "lib2_trimmed.fastq"],
    min_qual_score=20, min_read_len=30, verbose=False,
)
```

## 15) `dnakit.epigenome.methylpy_remove_clonal_bam()`

- **作用：** 用 Picard MarkDuplicates 从 BAM 文件中去除 PCR 重复（克隆）读段。
- **处理方法：** 检查输入存在、输出尚不存在，要求提供 `path_to_picard`，检查 `<path_to_picard>/picard.jar` 存在且 `java` 可用，然后调用 `methylpy.call_mc_se.remove_clonal_bam`，它运行 `java <java_options> -jar <path_to_picard>/picard.jar MarkDuplicates`（`ASSUME_SORTED=true`、`REMOVE_DUPLICATES=true`，宽松校验），写出重复度量文件 `metric` 并从中读回重复数（DNAKit 不返回该数值）。
- **API：** `dnakit.epigenome.methylpy_remove_clonal_bam(input_bam: 'str', output_bam: 'str', metric: 'str', *, is_pe: 'bool | None' = None, path_to_picard: 'str | None' = None, java_options: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`）、PATH 中的 `java`，以及包含 `picard.jar` 的文件夹（通过 `path_to_picard` 传入）。未给出该文件夹时抛出 `ConfigurationError`；文件夹中没有 `picard.jar` 时抛出 `BackendUnavailableError`。
- **输入：** `input_bam`：已存在的 BAM；`output_bam`：新路径；`metric`：重复度量文件；`path_to_picard`：`picard.jar` 所在文件夹（必填）；可选 `is_pe`（双端数据时设置）与 `java_options`（上游默认 `-Xmx20g`）。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（调用后存在的已声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_remove_clonal_bam
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

Path("sample.bam").touch()
for options in ({}, {"path_to_picard": "no_such_dir"}):
    try:
        methylpy_remove_clonal_bam("sample.bam", "dedup.bam", "metrics.txt", verbose=False, **options)
    except (ConfigurationError, BackendUnavailableError) as error:
        print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
ConfigurationError - path_to_picard must name the folder that contains picard.jar (methylpy reads path_to_picard + '/picard.jar').
BackendUnavailableError - picard.jar was not found at no_such_dir/picard.jar.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_remove_clonal_bam

result = methylpy_remove_clonal_bam(
    "s1_lib1.bam", "s1_lib1_no_clonal.bam", "s1_lib1_metrics.txt",
    path_to_picard="/opt/picard", is_pe=True, verbose=False,
)
print(result["outputs"])
```

## 16) `dnakit.epigenome.methylpy_run_rms_tests()`

- **作用：** 对准备好的按染色体文件运行 DMRfind 的均方根置换检验。
- **处理方法：** 检查内置的 `run_rms_tests.out` 二进制文件与输入文件存在，然后调用 `methylpy.DMRfind.run_rms_tests`，它在这些文件上启动该二进制文件并把检验结果写入 `output`。`DMRfind` 内部会调用这一步；只有已经持有 methylpy 中间文件时才需要直接使用它。
- **API：** `dnakit.epigenome.methylpy_run_rms_tests(files: 'Sequence[str]', output: 'str', samples: 'Sequence[str]', *, min_cov: 'int | None' = None, num_sims: 'int | str | None' = None, num_sig_tests: 'int | str | None' = None, seed: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`）；内置的 `run_rms_tests.out` 二进制文件需要 GSL 共享库 `libgsl.so.0` 和 `libgslcblas.so.0`（`apt install libgsl-dev` 或 `conda install -c conda-forge gsl`）。
- **输入：** `files`：已存在的输入文件列表；`output`：结果路径；`samples`：样本名。可选：`min_cov`、`num_sims` 与 `num_sig_tests`（整数或字符串）、`seed`、`keep_temp_files`。
- **输出：** 上述结果字典（`outputs` 为空列表）。输入文件不存在时抛出 `ConfigurationError`。
- **示例代码：**

```python
from dnakit.epigenome import methylpy_run_rms_tests
from dnakit.exceptions import ConfigurationError

try:
    methylpy_run_rms_tests(["missing_input.tsv"], "rms_out.tsv", ["a", "b"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, "-", error.message.split(":")[0])
```

- **示例输出：**

```text
ConfigurationError - files must be an existing file
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_run_rms_tests

result = methylpy_run_rms_tests(
    ["dmr_rms_input_chr1.tsv"], "dmr_rms_out.tsv", ["a", "b"], num_sims=1000, verbose=False
)
```

## 17) `dnakit.epigenome.methylpy_se_build_ref()`

- **作用：** 构建 C 转 T 与 G 转 A 的转换参考 FASTA 文件及其比对索引。
- **处理方法：** 检查输入与索引构建程序（`bowtie2-build`，`aligner="bowtie"` 时为 `bowtie-build`），然后调用 `methylpy.call_mc_se.build_ref`。上游按 `buffsize` 字节块读取 FASTA，写出 `<output>_f.fasta`（所有 C 替换为 T）和 `<output>_r.fasta`（所有 G 替换为 A），再为每个文件建立索引。双端函数同样使用它。
- **API：** `dnakit.epigenome.methylpy_se_build_ref(input_files: 'Sequence[str]', output: 'str', *, buffsize: 'int | None' = None, aligner: 'str | None' = None, path_to_aligner: 'str | None' = None, num_procs: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `bowtie2-build`（来自 bowtie2；`aligner="bowtie"` 时为 `bowtie-build`）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `input_files`：非空的已存在 FASTA 文件列表；`output`：两个转换参考序列的前缀。可选：`buffsize`、`aligner`（`bowtie2` 或 `bowtie`）、`path_to_aligner`（文件夹）、`num_procs`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_build_ref
from dnakit.exceptions import BackendUnavailableError

try:
    methylpy_se_build_ref(["ref.fa"], "genome", path_to_aligner="no_such_dir", verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - bowtie2-build was not found at no_such_dir/bowtie2-build.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_build_ref

result = methylpy_se_build_ref(["genome.fa"], "genome", num_procs=4, verbose=False)
# genome_f.fasta, genome_r.fasta and the bowtie2 indexes are written
```

## 18) `dnakit.epigenome.methylpy_se_call_methylated_sites()`

- **作用：** 从排序后的单端 BAM 文件检出甲基化胞嘧啶。
- **处理方法：** 检查 BAM 与参考序列存在，检查 `samtools`（`bgzip=True` 时还检查 `bgzip` 与 `tabix`），然后调用 `methylpy.call_mc_se.call_methylated_sites`。上游按染色体生成 pileup，用 `unmethylated_control` 估计未转化率，并以 `sig_cutoff` 调整后的 p 值做二项检验，每条染色体写出一个 `allc_<sample>_<chrom>.tsv`。
- **API：** `dnakit.epigenome.methylpy_se_call_methylated_sites(inputf: 'str', sample: 'str', reference_fasta: 'str', *, unmethylated_control: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_mpileup_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, min_mapq: 'int | None' = None, min_cov: 'int | None' = None, binom_test: 'bool | None' = None, path_to_samtools: 'str | None' = None, remove_chr_prefix: 'bool | None' = None, sort_mem: 'str | None' = None, add_snp_info: 'bool | None' = None, min_base_quality: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（`bgzip=True` 时还有 `bgzip` 与 `tabix`）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `inputf`：排序后的单端 BAM；`sample`：输出文件中使用的名称；`reference_fasta`：已存在的 FASTA。可选：`unmethylated_control`（`chrom:start-end` 区域或一个数值）、`sig_cutoff`、`num_procs`、`num_upstr_bases`、`num_downstr_bases`、`min_mapq`、`min_cov`、`min_base_quality`、`binom_test`、`generate_mpileup_file`、`compress_output`、`bgzip`、`remove_chr_prefix`、`add_snp_info`、`sort_mem`、`path_to_samtools`、`path_to_bgzip`、`path_to_tabix`、`keep_temp_files`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_call_methylated_sites
from dnakit.exceptions import BackendUnavailableError

Path("reads_se.bam").touch()
try:
    methylpy_se_call_methylated_sites(
        "reads_se.bam", "s1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_call_methylated_sites

result = methylpy_se_call_methylated_sites(
    "s1_processed_reads_no_clonal.bam", "s1", "genome.fa",
    unmethylated_control="chrC:", num_procs=4, verbose=False,
)
print(result["status"])  # allc_s1_<chrom>.tsv.gz files are written
```

## 19) `dnakit.epigenome.methylpy_se_call_methylated_sites_with_snp_info()`

- **作用：** 对单端 BAM 检出甲基化胞嘧啶，并在每个覆盖位点写出 SNP 信息。
- **处理方法：** 检查 BAM 与参考序列存在，检查 `samtools`（`bgzip=True` 时还检查 `bgzip` 与 `tabix`），然后调用 `methylpy.call_mc_se.call_methylated_sites_with_SNP_info`。它与 `methylpy_se_call_methylated_sites` 类似，但 allc 行有九列：常规的七列（最后一列固定为 1，不做二项检验），再加上上下文每个位置与预期转化不一致的碱基及其覆盖度。
- **API：** `dnakit.epigenome.methylpy_se_call_methylated_sites_with_snp_info(inputf: 'str', sample: 'str', reference_fasta: 'str', *, unmethylated_control: 'float | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_mpileup_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, min_mapq: 'int | None' = None, min_cov: 'int | None' = None, binom_test: 'bool | None' = None, path_to_samtools: 'str | None' = None, remove_chr_prefix: 'bool | None' = None, sort_mem: 'str | None' = None, add_snp_info: 'bool | None' = None, min_base_quality: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（`bgzip=True` 时还有 `bgzip` 与 `tabix`）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `inputf`：排序后的单端 BAM；`sample`：输出文件中使用的名称；`reference_fasta`：已存在的 FASTA。可选：`unmethylated_control`（此签名中类型为数值）、`sig_cutoff`、`num_procs`、`num_upstr_bases`、`num_downstr_bases`、`min_mapq`、`min_cov`、`min_base_quality`、`binom_test`、`generate_mpileup_file`、`compress_output`、`bgzip`、`remove_chr_prefix`、`add_snp_info`、`sort_mem`、`path_to_samtools`、`path_to_bgzip`、`path_to_tabix`、`keep_temp_files`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_call_methylated_sites_with_snp_info
from dnakit.exceptions import BackendUnavailableError

Path("reads_se.bam").touch()
try:
    methylpy_se_call_methylated_sites_with_snp_info(
        "reads_se.bam", "s1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_call_methylated_sites_with_snp_info

result = methylpy_se_call_methylated_sites_with_snp_info(
    "s1_processed_reads_no_clonal.bam", "s1", "genome.fa", num_procs=4, verbose=False
)
```

## 20) `dnakit.epigenome.methylpy_se_find_multi_mappers()`

- **作用：** 从单端比对中筛出质量足够的读段，并还原原始读段序列。
- **处理方法：** 调用 `methylpy.call_mc_se.find_multi_mappers`。它读取比对程序写出的 SAM 文件，跳过头部行、未比对读段、比对质量低于 `max(3, min_mapq)` 的读段以及补充比对，并利用读段名中 `!` 之后编码的位置还原亚硫酸氢盐转换前的序列（反向链读段先反转）。保留的行被轮流写入 `<output>_sorted_<i>`（`i < num_procs`）。除非 `keep_temp_files=True`，输入 SAM 文件会被删除。
- **API：** `dnakit.epigenome.methylpy_se_find_multi_mappers(inputf: 'str', output: 'str', *, num_procs: 'int | None' = None, min_mapq: 'int | None' = None, keep_temp_files: 'bool | None' = None, append: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 只需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），不运行外部程序。输入必须是读段名以 `!<编码位置>` 结尾的 SAM 文件，即 methylpy 比对步骤写出的格式。
- **输入：** `inputf`：已存在的 SAM 文件；`output`：结果文件前缀。可选：`num_procs`（输出文件数）、`min_mapq`、`keep_temp_files`、`append`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（此处恒为空列表，因为该函数没有声明单个输出路径：methylpy 按前缀命名并写出文件）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_se_find_multi_mappers

folder = Path(tempfile.mkdtemp())
sam = folder / "reads.sam"
sam.write_text(
    "@HD\tVN:1.0\n"
    "r1!\t0\tchr1\t10\t42\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
    "r2!\t16\tchr1\t50\t42\t8M\t*\t0\t0\tTTGGTTAA\tIIIIIIII\n"
    "r3!\t0\tchr1\t80\t1\t8M\t*\t0\t0\tCCGGCCGG\tIIIIIIII\n"
    "r4!\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAA\tIIIIIIII\n"
)
result = methylpy_se_find_multi_mappers(str(sam), str(folder / "out"), verbose=False)
rows = (folder / "out_sorted_0").read_text().splitlines()
print(result["status"], result["outputs"], sam.exists())
print([row.split("\t")[:2] for row in rows])
```

- **示例输出：**

```text
completed [] False
[['r1', '0'], ['r2', '16']]
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_find_multi_mappers

result = methylpy_se_find_multi_mappers("s1_forward.sam", "s1_forward", num_procs=4, verbose=False)
# s1_forward_sorted_0 ... s1_forward_sorted_3 are written; s1_forward.sam is deleted
```

## 21) `dnakit.epigenome.methylpy_se_merge_bam_files()`

- **作用：** 把多个单端 BAM 文件合并为一个，并修正头部。
- **处理方法：** 检查输入存在且 `samtools` 可用，然后调用 `methylpy.call_mc_se.merge_bam_files`，它在当前目录写出临时的 `header.sam`（第一个文件的头部加每个输入一行 `@RG`，随后删除），并运行 `samtools merge -r -h header.sam`。每个输入文件名都必须包含 `.bam`。
- **API：** `dnakit.epigenome.methylpy_se_merge_bam_files(input_files: 'Sequence[str]', output: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `input_files`：非空的已存在 BAM 文件列表；`output`：合并后的 BAM；`path_to_samtools`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（调用后存在的已声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_merge_bam_files
from dnakit.exceptions import BackendUnavailableError

Path("lib1.bam").touch()
Path("lib2.bam").touch()
try:
    methylpy_se_merge_bam_files(
        ["lib1.bam", "lib2.bam"], "merged.bam", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_merge_bam_files

result = methylpy_se_merge_bam_files(
    ["s1_lib1_no_clonal.bam", "s1_lib2_no_clonal.bam"], "s1_merged.bam", verbose=False
)
print(result["outputs"])
```

## 22) `dnakit.epigenome.methylpy_se_merge_sorted_multimap()`

- **作用：** 合并一个文库的排序后单端比对文件，只保留唯一比对的读段。
- **处理方法：** 检查文件与参考序列存在且 `samtools` 可用，然后调用 `methylpy.call_mc_se.merge_sorted_multimap`。当 `<reference_fasta>.fai` 不存在时，上游用 `samtools faidx` 为参考序列建索引，根据索引写出 SAM 头部，逐条合并排序后的文件并只保留在所有文件中仅出现一次的读段，最后转换为 `<prefix>_processed_reads.bam`。
- **API：** `dnakit.epigenome.methylpy_se_merge_sorted_multimap(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `current_library`：文库 ID；`files`：`methylpy_se_find_multi_mappers` 写出的 `_sorted_<i>` 文件；`prefix`：输出前缀；`reference_fasta`：已存在的 FASTA。`path_to_samtools` 可选。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_merge_sorted_multimap
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_se_merge_sorted_multimap(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_merge_sorted_multimap

result = methylpy_se_merge_sorted_multimap(
    "lib1", ["s1_lib1_sorted_0", "s1_lib1_sorted_1"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 23) `dnakit.epigenome.methylpy_se_merge_sorted_multimap_max_mapq()`

- **作用：** 合并一个文库的排序后单端比对文件，并为多重比对读段保留 MAPQ 最高的比对。
- **处理方法：** 检查文件与参考序列存在且 `samtools` 可用，然后调用 `methylpy.call_mc_se.merge_sorted_multimap_max_mapq`。当 `<reference_fasta>.fai` 不存在时，上游用 `samtools faidx` 为参考序列建索引，根据索引写出 SAM 头部，逐条合并排序后的文件并对出现多次的读段保留比对质量最高的比对，最后转换为 `<prefix>_processed_reads.bam`。
- **API：** `dnakit.epigenome.methylpy_se_merge_sorted_multimap_max_mapq(current_library: 'str', files: 'Sequence[str]', prefix: 'str', reference_fasta: 'str', *, path_to_samtools: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** 与 `methylpy_se_merge_sorted_multimap` 相同：`current_library`、`files`、`prefix`、已存在的 `reference_fasta`，以及可选的 `path_to_samtools`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_merge_sorted_multimap_max_mapq
from dnakit.exceptions import BackendUnavailableError

Path("lib1_sorted_0").touch()
try:
    methylpy_se_merge_sorted_multimap_max_mapq(
        "lib1", ["lib1_sorted_0"], "s1_lib1", "ref.fa", path_to_samtools="no_such_dir", verbose=False
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_merge_sorted_multimap_max_mapq

result = methylpy_se_merge_sorted_multimap_max_mapq(
    "lib1", ["s1_lib1_sorted_0"], "s1_lib1", "genome.fa", verbose=False
)
# s1_lib1_processed_reads.bam is written
```

## 24) `dnakit.epigenome.methylpy_se_run_alignment()`

- **作用：** 把一个文库的单端读段比对到两个转换参考序列上，并合并两条链的结果。
- **处理方法：** 检查文件与程序，然后调用 `methylpy.call_mc_se.run_alignment`。上游用所选比对程序（默认 `bowtie2`，也可用 `bowtie`）把读段分别比对到 C 转 T 参考序列和 G 转 A 参考序列，合并两条链并保留只在一处比对上的读段。
- **API：** `dnakit.epigenome.methylpy_se_run_alignment(current_library: 'str', library_read_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, aligner: 'str | None' = None, path_to_aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, num_procs: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools` 与 `bowtie2`（`aligner="bowtie"` 时为 `bowtie`）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `current_library`：文库 ID；`library_read_files`：已存在的 FASTQ 文件列表；`sample`；`forward_reference` 与 `reverse_reference`：`methylpy_se_build_ref` 构建的转换 FASTA（必须存在）；`reference_fasta`。可选：`path_to_output`、`aligner`、`aligner_options`、`path_to_aligner`、`path_to_samtools`、`merge_by_max_mapq`、`min_mapq`、`num_procs`、`keep_temp_files`、`sort_mem`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_run_alignment
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_se_run_alignment(
        "lib1", ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_run_alignment

result = methylpy_se_run_alignment(
    "lib1", ["lib1.fastq.gz"], "s1", "genome_f.fasta", "genome_r.fasta", "genome.fa",
    num_procs=4, verbose=False,
)
```

## 25) `dnakit.epigenome.methylpy_se_run_mapping()`

- **作用：** 对一个文库的单端亚硫酸氢盐读段完成修剪、比对与合并。
- **处理方法：** 检查文件与程序，然后调用 `methylpy.call_mc_se.run_mapping`。上游用所选比对程序（默认 `bowtie2`，也可用 `bowtie`）把读段分别比对到 C 转 T 参考序列和 G 转 A 参考序列，在可选的 cutadapt 接头修剪（`trim_reads`，默认开启）之后，合并两条链并保留只在一处比对上的读段。
- **API：** `dnakit.epigenome.methylpy_se_run_mapping(current_library: 'str', library_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, path_to_output: 'str | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, pbat: 'bool | None' = None, num_procs: 'int | None' = None, trim_reads: 'bool | None' = None, path_to_cutadapt: 'str | None' = None, adapter_seq: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, sort_mem: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`、`bowtie2`（或 `bowtie`）与 `cutadapt`（`trim_reads=False` 时不需要）（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** 与 `methylpy_se_run_alignment` 相同（用 `library_files` 代替 `library_read_files`；参考文件必须存在），另有 `pbat`、`trim_reads`、`path_to_cutadapt`、`adapter_seq`、`min_qual_score` 与 `min_read_len`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_run_mapping
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_se_run_mapping(
        "lib1", ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_run_mapping

result = methylpy_se_run_mapping(
    "lib1", ["lib1.fastq.gz"], "s1", "genome_f.fasta", "genome_r.fasta", "genome.fa",
    num_procs=4, verbose=False,
)
```

## 26) `dnakit.epigenome.methylpy_se_run_methylation_pipeline()`

- **作用：** 运行从 FASTQ 文件到 allc 表的完整单端亚硫酸氢盐流程。
- **处理方法：** 检查输入与程序，然后调用 `methylpy.call_mc_se.run_methylation_pipeline`。上游依次修剪接头，把每个文库比对到两个转换参考序列，合并文库，可选地用 Picard 去除克隆读段，并按染色体检出甲基化位点。这是串联本页其他单端函数的入口。
- **API：** `dnakit.epigenome.methylpy_se_run_methylation_pipeline(read_files: 'Sequence[str]', sample: 'str', forward_reference: 'str', reverse_reference: 'str', reference_fasta: 'str', *, libraries: 'Sequence[str] | None' = None, unmethylated_control: 'str | None' = None, path_to_output: 'str | None' = None, sig_cutoff: 'float | None' = None, num_procs: 'int | None' = None, sort_mem: 'str | None' = None, num_upstr_bases: 'int | None' = None, num_downstr_bases: 'int | None' = None, generate_allc_file: 'bool | None' = None, compress_output: 'bool | None' = None, bgzip: 'bool | None' = None, binom_test: 'bool | None' = None, min_cov: 'int | None' = None, trim_reads: 'bool | None' = None, aligner: 'str | None' = None, aligner_options: 'Sequence[str] | None' = None, merge_by_max_mapq: 'bool | None' = None, min_mapq: 'int | None' = None, remove_clonal: 'bool | None' = None, remove_chr_prefix: 'bool | None' = None, pbat: 'bool | None' = None, path_to_samtools: 'str | None' = None, path_to_aligner: 'str | None' = None, path_to_cutadapt: 'str | None' = None, path_to_picard: 'str | None' = None, path_to_bgzip: 'str | None' = None, path_to_tabix: 'str | None' = None, adapter_seq: 'str | None' = None, min_qual_score: 'int | None' = None, min_read_len: 'int | None' = None, keep_temp_files: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），以及 PATH 中的外部程序 `samtools`、`bowtie2`（或 `bowtie`）、`cutadapt`（`trim_reads=False` 时不需要）、`bgzip` 与 `tabix`（仅 `bgzip=True` 时）；`remove_clonal=True` 还需要 Picard 与 Java，它们不会被预先检查（或通过对应的 `path_to_...` 参数指定）。缺少时调用会抛出 `BackendUnavailableError`。
- **输入：** `read_files`：已存在的 FASTQ 文件列表；`sample`；`forward_reference` 与 `reverse_reference`（已存在的转换 FASTA）；`reference_fasta`。可选：`libraries`、`unmethylated_control`、`path_to_output`、`sig_cutoff`、`num_procs`、`min_cov`、`trim_reads`、`aligner`、`aligner_options`、`merge_by_max_mapq`、`min_mapq`、`remove_clonal`、`remove_chr_prefix`、`pbat`、`generate_allc_file`、`compress_output`、`bgzip`、`binom_test`、各 `path_to_...` 参数、`adapter_seq`、`min_qual_score`、`min_read_len` 与 `keep_temp_files`。
- **输出：** 字典，包含 `status`（`"completed"`）、`backend`（`"methylpy"`）、methylpy 的 `version`、`function`（实际运行的 methylpy 函数）、`outputs`（空列表：该函数没有声明输出路径）和 `arguments`（被转发的参数名，已排序）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.epigenome import methylpy_se_run_methylation_pipeline
from dnakit.exceptions import BackendUnavailableError

Path("ref_f.fasta").touch()
Path("ref_r.fasta").touch()
try:
    methylpy_se_run_methylation_pipeline(
        ["reads.fq"], "s1", "ref_f.fasta", "ref_r.fasta", "ref.fa",
        path_to_samtools="no_such_dir", verbose=False,
    )
except BackendUnavailableError as error:
    print(type(error).__name__, "-", error.message)
```

- **示例输出：**

```text
BackendUnavailableError - samtools was not found at no_such_dir/samtools.
```

安装依赖后的用法（不执行）：

```python
from dnakit.epigenome import methylpy_se_run_methylation_pipeline

result = methylpy_se_run_methylation_pipeline(
    ["lib1.fastq.gz"], "s1", "genome_f.fasta", "genome_r.fasta", "genome.fa",
    unmethylated_control="chrC:", num_procs=8, verbose=False,
)
```

## 27) `dnakit.epigenome.methylpy_filter_allc_files()`

- **作用：** 按胞嘧啶上下文、染色体和读段覆盖度过滤 allc 表。
- **处理方法：** 校验参数并统计输入行数，然后顺序调用 `methylpy.utilities.filter_allc_files`（`num_procs` 固定为 1）并捕获其输出。当一行的染色体在 `chroms` 中（若给出）、上下文匹配某个 `mc_type` 模式（`H` 与 `N` 为通配符）且 `min_cov <= 总读段数（第 6 列） <= max_cov` 时保留该行。`mc_type` 传入裸字符串会被拒绝，因为 methylpy 会按子串匹配并保留所有胞嘧啶。methylpy 的 SNP 错配过滤未提供，因为在 methylpy 1.4.7 中它们在 Python 3 上不起作用。已存在的输出不会被覆盖，失败时会删除部分输出。
- **API：** `dnakit.epigenome.methylpy_filter_allc_files(allc_files: 'Sequence[str | Path]', output_files: 'Sequence[str | Path]', *, mc_type: 'Sequence[str] | None' = None, chroms: 'Sequence[str] | None' = None, min_cov: 'int' = 0, max_cov: 'int | None' = None, compress_output: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 只需要 methylpy（额外依赖 `dnakit[bio-methylpy]`），不运行外部程序。
- **输入：** `allc_files`：1 至 32 个已存在的 allc 表（纯文本或 gzip）；`output_files`：每个输入对应一个新路径（`.gz` 名称需要 `compress_output=True`；启用时若缺少 `.gz` 则 methylpy 会追加）。`mc_type`：三字母上下文列表，如 `["CGN"]`；`chroms`：染色体名列表；`min_cov`（默认 0）与 `max_cov`（默认无上限），均为闭区间；`compress_output`（默认 False）。
- **输出：** 字典，包含 `status`、`operation`、`inputs`、`outputs`（每个文件的 `input`、`output`、`records_in`、`records_out`）、`filters`（实际使用的 `mc_type`、`chroms`、`min_cov`、`max_cov`）、`compressed` 和 `provenance`（methylpy 版本、使用的规则和未提供的错配过滤）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.epigenome import methylpy_filter_allc_files

output = Path(tempfile.mkdtemp()) / "cg.tsv"
result = methylpy_filter_allc_files(
    ["allc.tsv"], [output], mc_type=["CGN"], min_cov=5, verbose=False
)
row = result["outputs"][0]
print(row["records_in"], row["records_out"], result["filters"]["mc_type"])
for line in output.read_text().splitlines():
    print(line)
```

- **示例输出：**

```text
5 2 ['CGN']
chr1	100	+	CGG	10	12	1
chr2	50	+	CGT	0	8	0
```
