# GATK 命令行封装：拷贝数、读段计数、覆盖度与 pileup

对 `gatk` 命令行的带类型封装，涵盖体细胞与胚系拷贝数流程（区间分箱、注释与过滤、读段与等位基因计数、正常样本面板、去噪、分段、调用、绘图、胚系 CNV 调用）、等位基因特异性读段计数、可调用位点、覆盖深度、表达计数以及 Mutect2 的辅助表。 每个函数构造固定的 `gatk` 命令行（不经 shell），只写入新路径，并返回可 JSON 序列化的字典；GATK 采用 BSD-3-Clause 许可（Picard 为 MIT），不随包分发。函数位于 `dnakit.variants.gatk_p00_copynumber_cli_functions`，并从 `dnakit.variants` 导出；示例从模块导入，因此两种情况下都能运行。

## 示例数据

所有示例共用第一个示例写入当前目录的小数据集：双染色体参考序列（`toy.fa`，1,000 bp）及 `.fai` 与 `.dict`、含 38 条双端读段并带一个读组（样本 `S1`）的已索引 BAM 与 CRAM、带 tabix 索引的 15 个杂合位点 VCF、BED 与 interval list、单倍型图谱、refFlat 与 GFF3 注释、chain 文件以及带 UMI 标签的 BAM。其余示例把输出写入临时目录并传入 `verbose=False`。`gatk` 启动脚本从 PATH（或 `DNAKIT_GATK`）查找；绘制图表的示例还需要 PATH 中有 `Rscript`。

**示例代码：**

```python
import random
from pathlib import Path

import pysam

rng = random.Random(7)
chr1 = "".join(rng.choice("ACGT") for _ in range(600))
chr2 = "".join(rng.choice("ACGT") for _ in range(400))
contigs = {"chr1": chr1, "chr2": chr2}
with open("toy.fa", "w") as fasta:
    for name, sequence in contigs.items():
        fasta.write(f">{name}\n")
        fasta.writelines(sequence[i : i + 60] + "\n" for i in range(0, len(sequence), 60))
Path("toy.fa.fai").write_text("chr1\t600\t6\t60\t61\nchr2\t400\t622\t60\t61\n")
Path("toy.dict").write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:600\n@SQ\tSN:chr2\tLN:400\n")
Path("toy.bed").write_text("chr1\t100\t300\tt1\nchr2\t50\t200\tt2\n")
Path("toy.interval_list").write_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n@SQ\tSN:chr2\tLN:400\n"
    "chr1\t101\t300\t+\tt1\nchr2\t51\t200\t+\tt2\n"
)
header = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 600}, {"SN": "chr2", "LN": 400}],
    "RG": [{"ID": "rg1", "SM": "S1", "LB": "lib1", "PL": "ILLUMINA", "PU": "unit1"}],
}
complement = str.maketrans("ACGT", "TGCA")
reads = []
for contig, sequence in contigs.items():
    for start in range(0, len(sequence) - 120, 40):
        name = f"r{len(reads) // 2 + 1}"
        reads.append((contig, start, False, sequence[start : start + 50], name))
        mate = sequence[start + 50 : start + 100][::-1].translate(complement)
        reads.append((contig, start + 50, True, mate, name))
reads.sort(key=lambda item: (item[0], item[1]))
with pysam.AlignmentFile("toy.bam", "wb", header=header) as bam:
    for contig, start, reverse, sequence, name in reads:
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = sequence
        read.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + rng.randint(25, 40)) for _ in sequence))
        read.reference_id = 0 if contig == "chr1" else 1
        read.reference_start = start
        read.cigarstring = "50M"
        read.mapping_quality = 60
        read.flag = 147 if reverse else 99
        read.next_reference_id = read.reference_id
        read.next_reference_start = start + 50 if not reverse else start - 50
        read.template_length = -100 if reverse else 100
        read.set_tags([("RG", "rg1"), ("NM", 0)])
        bam.write(read)
pysam.index("toy.bam")
pysam.view("-C", "-T", "toy.fa", "-o", "toy.cram", "toy.bam", catch_stdout=False)
vcf = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=600>",
    "##contig=<ID=chr2,length=400>",
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
]
sites = []
for contig, sequence in contigs.items():
    for position in range(40, len(sequence) - 40, 60):
        ref = sequence[position - 1]
        alt = next(base for base in "ACGT" if base != ref)
        sites.append((contig, position, ref, alt))
        vcf.append(f"{contig}\t{position}\trs{position}{contig[-1]}\t{ref}\t{alt}\t50\tPASS\tAF=0.4\tGT\t0/1")
Path("toy.vcf").write_text("\n".join(vcf) + "\n")
pysam.tabix_index("toy.vcf", preset="vcf", force=True, keep_original=True)
Path("toy.haplotype_map").write_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:600\n@SQ\tSN:chr2\tLN:400\n"
    + "".join(f"{c}\t{p}\trs{p}{c[-1]}\t{r}\t{a}\t0.4\trs{p}{c[-1]}\n" for c, p, r, a in sites)
)
Path("toy.refFlat").write_text(
    "G1\tt1\tchr1\t+\t0\t600\t0\t600\t2\t0,300,\t200,600,\nG2\tt2\tchr2\t+\t0\t400\t0\t400\t1\t0,\t400,\n"
)
Path("toy.gff3").write_text(
    "##gff-version 3\nchr1\tsrc\tgene\t1\t600\t.\t+\t.\tID=G1;Name=G1\n"
    "chr1\tsrc\texon\t1\t200\t.\t+\t.\tParent=G1;Name=G1\nchr1\tsrc\texon\t301\t600\t.\t+\t.\tParent=G1;Name=G1\n"
    "chr2\tsrc\tgene\t1\t400\t.\t+\t.\tID=G2;Name=G2\nchr2\tsrc\texon\t1\t400\t.\t+\t.\tParent=G2;Name=G2\n"
)
Path("shift.chain").write_text("chain 1000 chr1 590 + 0 590 chr1 600 + 10 600 1\n590\n\n")
Path("old.interval_list").write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:590\nchr1\t101\t200\t+\tt1\n")
Path("old.haplotype_map").write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:590\nchr1\t30\trs1\tT\tA\t0.4\trs1\n")
with pysam.AlignmentFile("toy.bam") as source, pysam.AlignmentFile("umi.bam", "wb", template=source) as out:
    for index, read in enumerate(source):
        read.set_tag("RX", "ACGT"[index % 4] * 3 + "-" + "ACGT"[(index // 4) % 4] * 3)
        out.write(read)
pysam.index("umi.bam")
print("toy data written")
```

**示例输出：**

```text
toy data written
```

## 1) `dnakit.variants.gatk_preprocess_intervals()`

- **作用：** 为覆盖度收集准备等长分箱（或带填充的靶区），并写成 interval list。
- **处理方法：** 由带类型的参数构造 `gatk PreprocessIntervals` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk PreprocessIntervals --help` 校验。
- **API：** `dnakit.variants.gatk_preprocess_intervals(reference: 'str', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, bin_length: 'int | None' = None, padding: 'int | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`reference`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2908。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `interval_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_preprocess_intervals

work = Path(tempfile.mkdtemp())
result = gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
print(result["interval_count"])
```

- **示例输出：**

```text
10
```

## 2) `dnakit.variants.gatk_annotate_intervals()`

- **作用：** 为区间注释 GC 含量，并在提供轨道时注释可比对性与片段重复含量。
- **处理方法：** 由带类型的参数构造 `gatk AnnotateIntervals` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk AnnotateIntervals --help` 校验。
- **API：** `dnakit.variants.gatk_annotate_intervals(intervals: 'Sequence[str]', reference: 'str', output_path: 'str', *, mappability_track: 'str | None' = None, segmental_duplication_track: 'str | None' = None, feature_query_lookahead: 'int | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`intervals`, `reference`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2813。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `rows`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_annotate_intervals,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
bins = gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
result = gatk_annotate_intervals([work / "bins.interval_list"], "toy.fa", work / "annotated.tsv", verbose=False)
print(result["rows"][0])
```

- **示例输出：**

```text
{'CONTIG': 'chr1', 'START': 1, 'END': 100, 'GC_CONTENT': 0.42}
```

## 3) `dnakit.variants.gatk_collect_read_counts()`

- **作用：** 在指定区间收集读段计数（拷贝数分析的覆盖度输入），写成 TSV 或 HDF5 文件。
- **处理方法：** 由带类型的参数构造 `gatk CollectReadCounts` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectReadCounts --help` 校验。
- **API：** `dnakit.variants.gatk_collect_read_counts(input_paths: 'Sequence[str]', intervals: 'Sequence[str]', output_path: 'str', *, output_format: 'str | None' = None, reference: 'str | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `intervals`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2833。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `rows`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_collect_read_counts,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
result = gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
print([row["COUNT"] for row in result["rows"]])
```

- **示例输出：**

```text
[5, 4, 6, 4, 5, 0, 5, 4, 5, 0]
```

## 4) `dnakit.variants.gatk_collect_allelic_counts()`

- **作用：** 在 interval list 或 VCF 的位点收集参考与替代等位基因计数。
- **处理方法：** 由带类型的参数构造 `gatk CollectAllelicCounts` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectAllelicCounts --help` 校验。
- **API：** `dnakit.variants.gatk_collect_allelic_counts(input_paths: 'Sequence[str]', intervals: 'Sequence[str]', reference: 'str', output_path: 'str', *, minimum_base_quality: 'int | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `intervals`, `reference`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2831。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `rows`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_collect_allelic_counts

work = Path(tempfile.mkdtemp())
result = gatk_collect_allelic_counts(["toy.bam"], ["toy.vcf"], "toy.fa", work / "allelic.tsv", verbose=False)
print(len(result["rows"]), result["rows"][1])
```

- **示例输出：**

```text
15 {'CONTIG': 'chr1', 'POSITION': 100, 'REF_COUNT': 1, 'ALT_COUNT': 2, 'REF_NUCLEOTIDE': 'G', 'ALT_NUCLEOTIDE': 'C'}
```

## 5) `dnakit.variants.gatk_filter_intervals()`

- **作用：** 按注释（GC 含量、可比对性、片段重复）与计数统计过滤区间。
- **处理方法：** 由带类型的参数构造 `gatk FilterIntervals` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk FilterIntervals --help` 校验。
- **API：** `dnakit.variants.gatk_filter_intervals(intervals: 'Sequence[str]', output_path: 'str', *, annotated_intervals: 'str | None' = None, input_paths: 'Sequence[str] | None' = None, minimum_gc_content: 'float | None' = None, maximum_gc_content: 'float | None' = None, minimum_mappability: 'float | None' = None, maximum_mappability: 'float | None' = None, minimum_segmental_duplication_content: 'float | None' = None, maximum_segmental_duplication_content: 'float | None' = None, low_count_filter_count_threshold: 'int | None' = None, low_count_filter_percentage_of_samples: 'float | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`intervals`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2819。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `interval_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_annotate_intervals,
    gatk_filter_intervals,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_annotate_intervals([work / "bins.interval_list"], "toy.fa", work / "annotated.tsv", verbose=False)
result = gatk_filter_intervals([work / "bins.interval_list"], work / "kept.interval_list", annotated_intervals=work / "annotated.tsv", minimum_gc_content=0.45, maximum_gc_content=0.55, verbose=False)
print(result["interval_count"])
```

- **示例输出：**

```text
9
```

## 6) `dnakit.variants.gatk_create_read_count_panel_of_normals()`

- **作用：** 用正常样本的读段计数创建用于读段计数去噪的正常样本面板（PoN，HDF5）。
- **处理方法：** 由带类型的参数构造 `gatk CreateReadCountPanelOfNormals` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CreateReadCountPanelOfNormals --help` 校验。
- **API：** `dnakit.variants.gatk_create_read_count_panel_of_normals(input_paths: 'Sequence[str]', output_path: 'str', *, annotated_intervals: 'str | None' = None, number_of_eigensamples: 'int | None' = None, minimum_interval_median_percentile: 'float | None' = None, maximum_zeros_in_sample_percentage: 'float | None' = None, maximum_zeros_in_interval_percentage: 'float | None' = None, do_impute_zeros: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2816。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_create_read_count_panel_of_normals

work = Path(tempfile.mkdtemp())
try:
    gatk_create_read_count_panel_of_normals(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 7) `dnakit.variants.gatk_denoise_read_counts()`

- **作用：** 用正常样本面板（或注释区间）对读段计数去噪，得到标准化与去噪后的拷贝数比。
- **处理方法：** 由带类型的参数构造 `gatk DenoiseReadCounts` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk DenoiseReadCounts --help` 校验。
- **API：** `dnakit.variants.gatk_denoise_read_counts(input_path: 'str', standardized_copy_ratios: 'str', denoised_copy_ratios: 'str', *, count_panel_of_normals: 'str | None' = None, annotated_intervals: 'str | None' = None, number_of_eigensamples: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `standardized_copy_ratios`, `denoised_copy_ratios`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2817。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `rows`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_collect_read_counts,
    gatk_denoise_read_counts,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
result = gatk_denoise_read_counts(work / "counts.tsv", work / "standardized.tsv", work / "denoised.tsv", verbose=False)
print(result["rows"][:2])
```

- **示例输出：**

```text
[{'CONTIG': 'chr1', 'START': 1, 'END': 100, 'LOG2_COPY_RATIO': 0.160964}, {'CONTIG': 'chr1', 'START': 101, 'END': 200, 'LOG2_COPY_RATIO': -0.160964}]
```

## 8) `dnakit.variants.gatk_model_segments()`

- **作用：** 用去噪拷贝数比对拷贝数分段建模，并用等位基因计数对次要等位基因分数分段建模。
- **处理方法：** 由带类型的参数构造 `gatk ModelSegments` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk ModelSegments --help` 校验。
- **API：** `dnakit.variants.gatk_model_segments(output_dir: 'str', output_prefix: 'str', *, denoised_copy_ratios: 'str | None' = None, allelic_counts: 'str | None' = None, normal_allelic_counts: 'str | None' = None, minimum_total_allele_count_case: 'int | None' = None, minimum_total_allele_count_normal: 'int | None' = None, genotyping_homozygous_log_ratio_threshold: 'float | None' = None, genotyping_base_error_rate: 'float | None' = None, maximum_number_of_segments_per_chromosome: 'int | None' = None, window_sizes: 'Sequence[int] | None' = None, number_of_changepoints_penalty_factor: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`output_dir`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2823。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `segment_files`, `segments`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_collect_read_counts,
    gatk_denoise_read_counts,
    gatk_model_segments,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
gatk_denoise_read_counts(work / "counts.tsv", work / "standardized.tsv", work / "denoised.tsv", verbose=False)
result = gatk_model_segments(work / "segments", "S1", denoised_copy_ratios=work / "denoised.tsv", verbose=False)
print(len(result["segments"]), sorted(path.rsplit("/", 1)[1] for path in result["segment_files"]))
```

- **示例输出：**

```text
2 ['S1.af.igv.seg', 'S1.cr.igv.seg', 'S1.cr.seg', 'S1.modelBegin.af.param', 'S1.modelBegin.cr.param', 'S1.modelBegin.seg', 'S1.modelFinal.af.param', 'S1.modelFinal.cr.param', 'S1.modelFinal.seg']
```

## 9) `dnakit.variants.gatk_call_copy_ratio_segments()`

- **作用：** 把拷贝数比分段调用为扩增（+）、缺失（-）或拷贝数中性（0）。
- **处理方法：** 由带类型的参数构造 `gatk CallCopyRatioSegments` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CallCopyRatioSegments --help` 校验。
- **API：** `dnakit.variants.gatk_call_copy_ratio_segments(input_path: 'str', output_path: 'str', *, neutral_segment_copy_ratio_lower_bound: 'float | None' = None, neutral_segment_copy_ratio_upper_bound: 'float | None' = None, outlier_neutral_segment_copy_ratio_z_score_threshold: 'float | None' = None, calling_copy_ratio_z_score_threshold: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2814。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `segments`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_call_copy_ratio_segments,
    gatk_collect_read_counts,
    gatk_denoise_read_counts,
    gatk_model_segments,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_collect_read_counts(["toy.bam"], [work / "bins.interval_list"], work / "counts.tsv", output_format="TSV", verbose=False)
gatk_denoise_read_counts(work / "counts.tsv", work / "standardized.tsv", work / "denoised.tsv", verbose=False)
gatk_model_segments(work / "segments", "S1", denoised_copy_ratios=work / "denoised.tsv", verbose=False)
result = gatk_call_copy_ratio_segments(work / "segments" / "S1.cr.seg", work / "S1.called.seg", verbose=False)
print([row["CALL"] for row in result["segments"]])
```

- **示例输出：**

```text
[0, 0]
```

## 10) `dnakit.variants.gatk_combine_segment_breakpoints()`

- **作用：** 合并两个分段文件的断点，并用每个文件中选定的列注释所得区间。
- **处理方法：** 由带类型的参数构造 `gatk CombineSegmentBreakpoints` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CombineSegmentBreakpoints --help` 校验。
- **API：** `dnakit.variants.gatk_combine_segment_breakpoints(segments: 'Sequence[str]', columns_of_interest: 'Sequence[str]', output_path: 'str', *, labels: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`segments`, `columns_of_interest`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2815。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `segments`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_combine_segment_breakpoints

work = Path(tempfile.mkdtemp())
try:
    gatk_combine_segment_breakpoints(segments=["missing.dat"], columns_of_interest=["x"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] segments must hold exactly two files.
```

## 11) `dnakit.variants.gatk_merge_annotated_regions()`

- **作用：** 合并带注释的基因组区域。不给 annotations_to_match 时合并相接或重叠的区域；给出时合并距离不超过 max_merge_distance 且所列注释值完全相同的区域。
- **处理方法：** 由带类型的参数构造 `gatk MergeAnnotatedRegions / MergeAnnotatedRegionsByAnnotation` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk MergeAnnotatedRegions / MergeAnnotatedRegionsByAnnotation --help` 校验。
- **API：** `dnakit.variants.gatk_merge_annotated_regions(segments: 'str', output_path: 'str', *, reference: 'str | None' = None, annotations_to_match: 'Sequence[str] | None' = None, max_merge_distance: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`segments`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2821, G2822。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `regions`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import (
    gatk_annotate_intervals,
    gatk_merge_annotated_regions,
    gatk_preprocess_intervals,
)

work = Path(tempfile.mkdtemp())
gatk_preprocess_intervals("toy.fa", work / "bins.interval_list", bin_length=100, padding=0, verbose=False)
gatk_annotate_intervals([work / "bins.interval_list"], "toy.fa", work / "annotated.tsv", verbose=False)
(work / "annotated.seg").write_bytes((work / "annotated.tsv").read_bytes())
plain = gatk_merge_annotated_regions(work / "annotated.seg", work / "plain.seg", reference="toy.fa", verbose=False)
by_annotation = gatk_merge_annotated_regions(work / "annotated.seg", work / "by_gc.seg", reference="toy.fa", annotations_to_match=["GC_CONTENT"], verbose=False)
print(len(plain["regions"]), by_annotation["tool"])
```

- **示例输出：**

```text
10 MergeAnnotatedRegionsByAnnotation
```

## 12) `dnakit.variants.gatk_plot_denoised_copy_ratios()`

- **作用：** 绘制标准化与去噪拷贝数比的图，并给出中位数绝对偏差。
- **处理方法：** 由带类型的参数构造 `gatk PlotDenoisedCopyRatios` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk PlotDenoisedCopyRatios --help` 校验。 The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API：** `dnakit.variants.gatk_plot_denoised_copy_ratios(standardized_copy_ratios: 'str', denoised_copy_ratios: 'str', sequence_dictionary: 'str', output_dir: 'str', output_prefix: 'str', *, minimum_contig_length: 'int | None' = None, maximum_copy_ratio: 'float | None' = None, point_size_copy_ratio: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`standardized_copy_ratios`, `denoised_copy_ratios`, `sequence_dictionary`, `output_dir`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2824。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `plot_files`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_plot_denoised_copy_ratios

work = Path(tempfile.mkdtemp())
try:
    gatk_plot_denoised_copy_ratios(standardized_copy_ratios="missing.dat", denoised_copy_ratios="missing.dat", sequence_dictionary="missing.dat", output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] standardized_copy_ratios must be an existing file
```

## 13) `dnakit.variants.gatk_plot_modeled_segments()`

- **作用：** 绘制去噪并分段后的拷贝数比与次要等位基因分数估计图。
- **处理方法：** 由带类型的参数构造 `gatk PlotModeledSegments` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk PlotModeledSegments --help` 校验。 The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API：** `dnakit.variants.gatk_plot_modeled_segments(segments: 'str', sequence_dictionary: 'str', output_dir: 'str', output_prefix: 'str', *, denoised_copy_ratios: 'str | None' = None, allelic_counts: 'str | None' = None, minimum_contig_length: 'int | None' = None, maximum_copy_ratio: 'float | None' = None, point_size_copy_ratio: 'float | None' = None, point_size_allele_fraction: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`segments`, `sequence_dictionary`, `output_dir`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2825。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `plot_files`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_plot_modeled_segments

work = Path(tempfile.mkdtemp())
try:
    gatk_plot_modeled_segments(segments="missing.dat", sequence_dictionary="missing.dat", output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] segments must be an existing file
```

## 14) `dnakit.variants.gatk_tag_germline_events()`

- **作用：** 给已调用的肿瘤分段文件中与匹配正常样本的胚系事件重叠的分段打上标签。
- **处理方法：** 由带类型的参数构造 `gatk TagGermlineEvents` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk TagGermlineEvents --help` 校验。
- **API：** `dnakit.variants.gatk_tag_germline_events(segments: 'str', called_matched_normal_seg_file: 'str', reference: 'str', output_path: 'str', *, endpoint_padding: 'int | None' = None, reciprocal_threshold: 'float | None' = None, input_call_header: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`segments`, `called_matched_normal_seg_file`, `reference`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2827。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `segments`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_tag_germline_events

work = Path(tempfile.mkdtemp())
try:
    gatk_tag_germline_events(segments="missing.dat", called_matched_normal_seg_file="missing.dat", reference="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] segments must be an existing file
```

## 15) `dnakit.variants.gatk_determine_germline_contig_ploidy()`

- **作用：** 依据读段计数确定胚系样本的基线染色体倍性（队列模式建立倍性模型，个案模式应用已有模型）。
- **处理方法：** 由带类型的参数构造 `gatk DetermineGermlineContigPloidy` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk DetermineGermlineContigPloidy --help` 校验。 This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the GATK conda environment (gatkcondaenv.yml of the release zip).
- **API：** `dnakit.variants.gatk_determine_germline_contig_ploidy(input_paths: 'Sequence[str]', output_dir: 'str', output_prefix: 'str', *, contig_ploidy_priors: 'str | None' = None, model: 'str | None' = None, intervals: 'Sequence[str] | None' = None, mapping_error_rate: 'float | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_dir`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2818。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `calls_directories`, `model_directories`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_determine_germline_contig_ploidy

work = Path(tempfile.mkdtemp())
try:
    gatk_determine_germline_contig_ploidy(input_paths=["missing.dat"], output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 16) `dnakit.variants.gatk_germline_cnv_caller()`

- **作用：** 依据读段计数与染色体倍性结果调用胚系拷贝数变异（COHORT 模式学习去噪模型，CASE 模式应用已有模型）。
- **处理方法：** 由带类型的参数构造 `gatk GermlineCNVCaller` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GermlineCNVCaller --help` 校验。 This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the GATK conda environment (gatkcondaenv.yml of the release zip).
- **API：** `dnakit.variants.gatk_germline_cnv_caller(run_mode: 'str', input_paths: 'Sequence[str]', contig_ploidy_calls: 'str', output_dir: 'str', output_prefix: 'str', *, model: 'str | None' = None, annotated_intervals: 'str | None' = None, intervals: 'Sequence[str] | None' = None, max_copy_number: 'int | None' = None, p_alt: 'float | None' = None, p_active: 'float | None' = None, interval_merging_rule: 'str | None' = 'OVERLAPPING_ONLY', options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`run_mode`, `input_paths`, `contig_ploidy_calls`, `output_dir`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2820。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `calls_directories`, `model_directories`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_germline_cnv_caller

work = Path(tempfile.mkdtemp())
try:
    gatk_germline_cnv_caller(run_mode="COHORT", input_paths=["missing.dat"], contig_ploidy_calls="missing_dir", output_dir=work / "d", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 17) `dnakit.variants.gatk_postprocess_germline_cnv_calls()`

- **作用：** 对单个样本的 GermlineCNVCaller 输出后处理，生成区间 VCF、分段 VCF 与去噪拷贝数比。
- **处理方法：** 由带类型的参数构造 `gatk PostprocessGermlineCNVCalls` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk PostprocessGermlineCNVCalls --help` 校验。 This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the GATK conda environment (gatkcondaenv.yml of the release zip).
- **API：** `dnakit.variants.gatk_postprocess_germline_cnv_calls(calls_shard_paths: 'Sequence[str]', model_shard_paths: 'Sequence[str]', contig_ploidy_calls: 'str', output_genotyped_intervals: 'str', output_genotyped_segments: 'str', output_denoised_copy_ratios: 'str', *, sample_index: 'int | None' = None, allosomal_contigs: 'Sequence[str] | None' = None, autosomal_ref_copy_number: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`calls_shard_paths`, `model_shard_paths`, `contig_ploidy_calls`, `output_genotyped_intervals`, `output_genotyped_segments`, `output_denoised_copy_ratios`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2826。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `segment_record_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_postprocess_germline_cnv_calls

work = Path(tempfile.mkdtemp())
try:
    gatk_postprocess_germline_cnv_calls(calls_shard_paths=["missing_dir"], model_shard_paths=["missing_dir"], contig_ploidy_calls="missing_dir", output_genotyped_intervals=work / "o", output_genotyped_segments=work / "o", output_denoised_copy_ratios=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] calls_shard_paths must be an existing directory
```

## 18) `dnakit.variants.gatk_ase_read_counter()`

- **作用：** 在杂合位点生成过滤后的碱基计数表，用于等位基因特异性表达分析。
- **处理方法：** 由带类型的参数构造 `gatk ASEReadCounter` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk ASEReadCounter --help` 校验。
- **API：** `dnakit.variants.gatk_ase_read_counter(input_paths: 'Sequence[str]', variant_paths: 'Sequence[str]', output_path: 'str', *, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, output_format: 'str | None' = None, min_base_quality: 'int | None' = None, min_mapping_quality: 'int | None' = None, min_depth_of_non_filtered_base: 'int | None' = None, max_depth_per_sample: 'int | None' = None, count_overlap_reads_handling: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `variant_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2828。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `rows`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_ase_read_counter

work = Path(tempfile.mkdtemp())
result = gatk_ase_read_counter(["toy.bam"], ["toy.vcf.gz"], work / "ase.tsv", reference="toy.fa", verbose=False)
print(result["rows"][2]["position"], result["rows"][2]["refCount"], result["rows"][2]["altCount"])
```

- **示例输出：**

```text
160 1 1
```

## 19) `dnakit.variants.gatk_analyze_saturation_mutagenesis()`

- **作用：** 处理 MITESeq 或其他饱和突变实验的读段，得到各变异的密码子与氨基酸改变计数。
- **处理方法：** 由带类型的参数构造 `gatk AnalyzeSaturationMutagenesis` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk AnalyzeSaturationMutagenesis --help` 校验。
- **API：** `dnakit.variants.gatk_analyze_saturation_mutagenesis(input_paths: 'Sequence[str]', reference: 'str', orf: 'str', output_file_prefix: 'str', *, min_q: 'int | None' = None, min_mapq: 'int | None' = None, min_length: 'int | None' = None, min_flanking_length: 'int | None' = None, min_variant_obs: 'int | None' = None, paired_mode: 'bool | None' = None, codon_translation: 'str | None' = None, write_rejected_reads: 'bool | None' = None, find_large_deletions: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `reference`, `orf`, `output_file_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2829。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `output_files`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_analyze_saturation_mutagenesis

work = Path(tempfile.mkdtemp())
try:
    gatk_analyze_saturation_mutagenesis(input_paths=["missing.dat"], reference="missing.dat", orf="x", output_file_prefix=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 20) `dnakit.variants.gatk_callable_loci()`

- **作用：** 确定每个位点的可调用状态（CALLABLE、LOW_COVERAGE、NO_COVERAGE、POOR_MAPPING_QUALITY、EXCESSIVE_COVERAGE、REF_N）。
- **处理方法：** 由带类型的参数构造 `gatk CallableLoci` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CallableLoci --help` 校验。
- **API：** `dnakit.variants.gatk_callable_loci(input_paths: 'Sequence[str]', reference: 'str', output_path: 'str', summary_path: 'str', *, intervals: 'Sequence[str] | None' = None, min_depth: 'int | None' = None, max_depth: 'int | None' = None, min_mapping_quality: 'int | None' = None, min_base_quality: 'int | None' = None, min_depth_for_low_mapq: 'int | None' = None, max_low_mapq: 'int | None' = None, max_fraction_of_reads_with_low_mapq: 'float | None' = None, output_format: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `reference`, `output_path`, `summary_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2830。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_callable_loci

work = Path(tempfile.mkdtemp())
result = gatk_callable_loci(["toy.bam"], "toy.fa", work / "callable.bed", work / "callable.txt", min_depth=2, verbose=False)
print(result["summary"])
```

- **示例输出：**

```text
{'REF_N': 0, 'CALLABLE': 720, 'NO_COVERAGE': 120, 'LOW_COVERAGE': 160, 'EXCESSIVE_COVERAGE': 0, 'POOR_MAPPING_QUALITY': 0}
```

## 21) `dnakit.variants.gatk_collect_f1r2_counts()`

- **作用：** 为 Mutect2 方向偏倚混合模型过滤器收集各序列上下文的 F1R2 读段计数。
- **处理方法：** 由带类型的参数构造 `gatk CollectF1R2Counts` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectF1R2Counts --help` 校验。
- **API：** `dnakit.variants.gatk_collect_f1r2_counts(input_paths: 'Sequence[str]', reference: 'str', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, f1r2_min_bq: 'int | None' = None, f1r2_median_mq: 'int | None' = None, f1r2_max_depth: 'int | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `reference`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2832。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_collect_f1r2_counts

work = Path(tempfile.mkdtemp())
result = gatk_collect_f1r2_counts(["toy.bam"], "toy.fa", work / "f1r2.tar.gz", verbose=False)
print(result["status"], (work / "f1r2.tar.gz").stat().st_size > 0)
```

- **示例输出：**

```text
completed True
```

## 22) `dnakit.variants.gatk_count_reads_or_bases()`

- **作用：** 统计 SAM/BAM/CRAM 文件在 GATK 默认读段过滤之后的读段数（unit='reads'）或碱基数（unit='bases'）。
- **处理方法：** 由带类型的参数构造 `gatk CountReads / CountBases` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CountReads / CountBases --help` 校验。
- **API：** `dnakit.variants.gatk_count_reads_or_bases(input_paths: 'Sequence[str]', output_path: 'str', *, unit: 'str | None' = None, intervals: 'Sequence[str] | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2834, G2835。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `count`, `unit`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_count_reads_or_bases

work = Path(tempfile.mkdtemp())
reads = gatk_count_reads_or_bases(["toy.bam"], work / "reads.txt", verbose=False)
bases = gatk_count_reads_or_bases(["toy.bam"], work / "bases.txt", unit="bases", verbose=False)
print(reads["count"], bases["count"])
```

- **示例输出：**

```text
38 1900
```

## 23) `dnakit.variants.gatk_depth_of_coverage()`

- **作用：** 为读段数据生成覆盖度汇总：每位点、每区间、每基因和每样本的深度，以及覆盖度直方图。
- **处理方法：** 由带类型的参数构造 `gatk DepthOfCoverage` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk DepthOfCoverage --help` 校验。
- **API：** `dnakit.variants.gatk_depth_of_coverage(input_paths: 'Sequence[str]', reference: 'str', intervals: 'Sequence[str]', output_prefix: 'str', *, partition_types: 'Sequence[str] | None' = None, summary_coverage_thresholds: 'Sequence[int] | None' = None, min_base_quality: 'int | None' = None, count_type: 'str | None' = None, output_format: 'str | None' = None, calculate_coverage_over_genes: 'str | None' = None, omit_depth_output_at_each_base: 'bool | None' = None, omit_interval_statistics: 'bool | None' = None, omit_locus_table: 'bool | None' = None, omit_per_sample_statistics: 'bool | None' = None, print_base_counts: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `reference`, `intervals`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2836。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `output_files`, `sample_summary`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_depth_of_coverage

work = Path(tempfile.mkdtemp())
result = gatk_depth_of_coverage(["toy.bam"], "toy.fa", ["toy.interval_list"], work / "coverage", verbose=False)
row = result["sample_summary"][0]
print(row["sample_id"], row["total"], row["mean"])
```

- **示例输出：**

```text
S1 860 2.46
```

## 24) `dnakit.variants.gatk_gene_expression_evaluation()`

- **作用：** 依据比对到基因组的 RNA-seq 读段评估基因表达：按 GFF 中的基因、外显子等特征统计读段数。
- **处理方法：** 由带类型的参数构造 `gatk GeneExpressionEvaluation` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GeneExpressionEvaluation --help` 校验。
- **API：** `dnakit.variants.gatk_gene_expression_evaluation(input_paths: 'Sequence[str]', gff_file: 'str', output_path: 'str', *, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, grouping_types: 'Sequence[str] | None' = None, overlap_types: 'Sequence[str] | None' = None, feature_label_key: 'str | None' = None, multi_map_method: 'str | None' = None, multi_overlap_method: 'str | None' = None, read_strands: 'str | None' = None, unspliced: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `gff_file`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2837。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `rows`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_gene_expression_evaluation

work = Path(tempfile.mkdtemp())
try:
    gatk_gene_expression_evaluation(input_paths=["missing.dat"], gff_file="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 25) `dnakit.variants.gatk_get_normal_artifact_data()`

- **作用：** 收集用于训练 Mutect2 正常样本伪迹过滤器的数据（正常样本在候选位点的等位基因计数）。
- **处理方法：** 由带类型的参数构造 `gatk GetNormalArtifactData` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GetNormalArtifactData --help` 校验。
- **API：** `dnakit.variants.gatk_get_normal_artifact_data(input_paths: 'Sequence[str]', normal_samples: 'Sequence[str]', reference: 'str', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, error_prob: 'float | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `normal_samples`, `reference`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2838。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_get_normal_artifact_data

work = Path(tempfile.mkdtemp())
result = gatk_get_normal_artifact_data(["toy.bam"], ["S1"], "toy.fa", work / "normal_artifact.table", verbose=False)
print(result["status"], (work / "normal_artifact.table").stat().st_size > 0)
```

- **示例输出：**

```text
completed True
```

## 26) `dnakit.variants.gatk_get_pileup_summaries()`

- **作用：** 汇总 pileup 指标（常见双等位位点的 ref、alt 与其他 alt 计数），用于推断样本间污染。
- **处理方法：** 由带类型的参数构造 `gatk GetPileupSummaries` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GetPileupSummaries --help` 校验。
- **API：** `dnakit.variants.gatk_get_pileup_summaries(input_paths: 'Sequence[str]', variant: 'str', intervals: 'Sequence[str]', output_path: 'str', *, minimum_population_allele_frequency: 'float | None' = None, maximum_population_allele_frequency: 'float | None' = None, max_depth_per_sample: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `variant`, `intervals`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2839。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `rows`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_get_pileup_summaries

work = Path(tempfile.mkdtemp())
result = gatk_get_pileup_summaries(["toy.bam"], "toy.vcf.gz", ["toy.vcf.gz"], work / "pileup.table", maximum_population_allele_frequency=0.9, verbose=False)
print(len(result["rows"]), result["rows"][2])
```

- **示例输出：**

```text
15 {'contig': 'chr1', 'position': 160, 'ref_count': 1, 'alt_count': 1, 'other_alt_count': 0, 'allele_frequency': 0.4}
```

## 27) `dnakit.variants.gatk_local_assembler()`

- **作用：** 对小区域做局部组装以发现结构变异（GATK BETA 工具），输出支架序列或组装图。
- **处理方法：** 由带类型的参数构造 `gatk LocalAssembler` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk LocalAssembler --help` 校验。 GATK marks this tool BETA (work in progress).
- **API：** `dnakit.variants.gatk_local_assembler(input_paths: 'Sequence[str]', intervals: 'Sequence[str]', assembly_name: 'str', *, fasta_file: 'str | None' = None, gfa_file: 'str | None' = None, reference: 'str | None' = None, min_sv_size: 'int | None' = None, no_scaffolding: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `intervals`, `assembly_name`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2840。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_local_assembler

work = Path(tempfile.mkdtemp())
result = gatk_local_assembler(["toy.bam"], ["chr1:1-600"], "demo", fasta_file=work / "scaffolds.fa", gfa_file=work / "graph.gfa", reference="toy.fa", verbose=False)
print(result["status"], sorted(result["outputs"]))
```

- **示例输出：**

```text
completed ['fasta_file', 'gfa_file']
```

## 28) `dnakit.variants.gatk_analyze_covariates()`

- **作用：** 评估并比较碱基质量分数重校准（BQSR）表并绘图。
- **处理方法：** 由带类型的参数构造 `gatk AnalyzeCovariates` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk AnalyzeCovariates --help` 校验。 The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API：** `dnakit.variants.gatk_analyze_covariates(*, before_report_file: 'str | None' = None, after_report_file: 'str | None' = None, bqsr_recal_file: 'str | None' = None, plots_report_file: 'str | None' = None, intermediate_csv_file: 'str | None' = None, ignore_last_modification_times: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2841。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.variants.gatk_p00_copynumber_cli_functions import gatk_analyze_covariates

work = Path(tempfile.mkdtemp())
try:
    gatk_analyze_covariates(verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] Give at least one of before_report_file, after_report_file or bqsr_recal_file.
```
