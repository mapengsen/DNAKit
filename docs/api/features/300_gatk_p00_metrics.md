# GATK 命令行封装：Illumina 运行、读段指标、指纹与文件检查

对 `gatk` 命令行（GATK 4.6.2.0，内含 Picard 3.4.0）的带类型封装，涵盖 Illumina 运行处理（碱基识别指标、条形码提取、转换为 FASTQ 或未比对 BAM、接头标记）、比对与测序指标收集（比对汇总、杂交捕获、WGS、RNA-seq、RRBS、伪迹、UMI、变异检出）、指纹分析，以及 SAM/BAM/CRAM 文件的检查、查看与比较。 每个函数构造固定的 `gatk` 命令行（不经 shell），只写入新路径，并返回可 JSON 序列化的字典；GATK 采用 BSD-3-Clause 许可（Picard 为 MIT），不随包分发。函数位于 `dnakit.sequencing_qc.gatk_p00_metrics_cli_functions`，并从 `dnakit.sequencing_qc` 导出；示例从模块导入，因此两种情况下都能运行。

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

## 1) `dnakit.sequencing_qc.gatk_collect_illumina_basecalling_metrics()`

- **作用：** 为测序运行的一条泳道收集 Illumina 碱基识别指标（簇数、PF 簇数、各条形码计数）。
- **处理方法：** 由带类型的参数构造 `gatk CollectIlluminaBasecallingMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectIlluminaBasecallingMetrics --help` 校验。 Needs the BaseCalls directory of an Illumina run (BCL/CBCL, filter, locs or clocs files).
- **API：** `dnakit.sequencing_qc.gatk_collect_illumina_basecalling_metrics(basecalls_dir: 'str', lane: 'int', read_structure: 'str', *, output_path: 'str | None' = None, barcodes_dir: 'str | None' = None, barcode_file: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`basecalls_dir`, `lane`, `read_structure`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2807。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_illumina_basecalling_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_illumina_basecalling_metrics(basecalls_dir="missing_dir", lane=1, read_structure="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] basecalls_dir must be an existing directory
```

## 2) `dnakit.sequencing_qc.gatk_collect_illumina_lane_metrics()`

- **作用：** 从运行目录收集 Illumina 泳道与相位指标（簇密度、PF 簇、phasing 与 prephasing）。
- **处理方法：** 由带类型的参数构造 `gatk CollectIlluminaLaneMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectIlluminaLaneMetrics --help` 校验。 Needs the InterOp files of an Illumina run.
- **API：** `dnakit.sequencing_qc.gatk_collect_illumina_lane_metrics(run_directory: 'str', output_directory: 'str', output_prefix: 'str', *, read_structure: 'str | None' = None, file_extension: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`run_directory`, `output_directory`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2808。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`metrics_files`, `outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_illumina_lane_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_illumina_lane_metrics(run_directory="missing_dir", output_directory="missing_dir", output_prefix="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] run_directory must be an existing directory
```

## 3) `dnakit.sequencing_qc.gatk_extract_illumina_barcodes()`

- **作用：** 确定 Illumina 泳道中每条读段的条形码，写出 *_barcode.txt 文件与各条形码指标。
- **处理方法：** 由带类型的参数构造 `gatk ExtractIlluminaBarcodes` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk ExtractIlluminaBarcodes --help` 校验。 Needs the BaseCalls directory of an Illumina run.
- **API：** `dnakit.sequencing_qc.gatk_extract_illumina_barcodes(basecalls_dir: 'str', lanes: 'Sequence[int]', read_structure: 'str', *, barcodes: 'Sequence[str] | None' = None, barcode_file: 'str | None' = None, metrics_file: 'str | None' = None, output_dir: 'str | None' = None, max_mismatches: 'int | None' = None, min_mismatch_delta: 'int | None' = None, max_no_calls: 'int | None' = None, minimum_base_quality: 'int | None' = None, minimum_quality: 'int | None' = None, distance_mode: 'str | None' = None, compress_outputs: 'bool | None' = None, num_processors: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`basecalls_dir`, `lanes`, `read_structure`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2809。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_extract_illumina_barcodes

work = Path(tempfile.mkdtemp())
try:
    gatk_extract_illumina_barcodes(basecalls_dir="missing_dir", lanes=[1], read_structure="x", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] Give exactly one of barcodes or barcode_file.
```

## 4) `dnakit.sequencing_qc.gatk_illumina_basecalls_to_reads()`

- **作用：** 把 Illumina 碱基识别数据转换为 FASTQ（output_format='fastq'）或未比对的 SAM/BAM/CRAM（output_format='sam'）。
- **处理方法：** 由带类型的参数构造 `gatk IlluminaBasecallsToFastq / IlluminaBasecallsToSam` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk IlluminaBasecallsToFastq / IlluminaBasecallsToSam --help` 校验。 Needs the BaseCalls directory of an Illumina run.
- **API：** `dnakit.sequencing_qc.gatk_illumina_basecalls_to_reads(basecalls_dir: 'str', lanes: 'Sequence[int]', read_structure: 'str', run_barcode: 'str', parameters_file: 'str', output_format: 'str', *, output_prefix: 'str | None' = None, sequencing_center: 'str | None' = None, metrics_file: 'str | None' = None, barcodes_dir: 'str | None' = None, include_non_pf_reads: 'bool | None' = None, max_mismatches: 'int | None' = None, min_mismatch_delta: 'int | None' = None, max_no_calls: 'int | None' = None, minimum_quality: 'int | None' = None, first_tile: 'int | None' = None, tile_limit: 'int | None' = None, num_processors: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`basecalls_dir`, `lanes`, `read_structure`, `run_barcode`, `parameters_file`, `output_format`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2810, G2811。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `output_format`, `metrics`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_illumina_basecalls_to_reads

work = Path(tempfile.mkdtemp())
try:
    gatk_illumina_basecalls_to_reads(basecalls_dir="missing_dir", lanes=[1], read_structure="x", run_barcode="x", parameters_file="missing.dat", output_format="fastq", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] basecalls_dir must be an existing directory
```

## 5) `dnakit.sequencing_qc.gatk_mark_illumina_adapters()`

- **作用：** 在读段中查找 Illumina 接头序列，重写文件并添加 XT 接头修剪标签，同时输出接头剪切直方图。
- **处理方法：** 由带类型的参数构造 `gatk MarkIlluminaAdapters` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk MarkIlluminaAdapters --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_mark_illumina_adapters(input_path: 'str', metrics_file: 'str', *, output_path: 'str | None' = None, adapters: 'Sequence[str] | None' = None, five_prime_adapter: 'str | None' = None, three_prime_adapter: 'str | None' = None, max_error_rate_se: 'float | None' = None, max_error_rate_pe: 'float | None' = None, min_match_bases_se: 'int | None' = None, min_match_bases_pe: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `metrics_file`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2812。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_mark_illumina_adapters

work = Path(tempfile.mkdtemp())
try:
    gatk_mark_illumina_adapters(input_path="missing.dat", metrics_file=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 6) `dnakit.sequencing_qc.gatk_collect_alignment_summary_metrics()`

- **作用：** 从 SAM/BAM/CRAM 文件生成比对指标汇总（读段数、PF 读段、比对碱基、错配与插入缺失率、嵌合体、接头）。
- **处理方法：** 由带类型的参数构造 `gatk CollectAlignmentSummaryMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectAlignmentSummaryMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_alignment_summary_metrics(input_path: 'str', output_path: 'str', *, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, assume_sorted: 'bool | None' = None, is_bisulfite_sequenced: 'bool | None' = None, max_insert_size: 'int | None' = None, expected_pair_orientations: 'Sequence[str] | None' = None, adapter_sequences: 'Sequence[str] | None' = None, histogram_file: 'str | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2850。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_alignment_summary_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_alignment_summary_metrics("toy.bam", work / "summary.txt", reference="toy.fa", verbose=False)
row = next(r for r in result["metrics"] if r["CATEGORY"] == "PAIR")
print(result["status"], row["TOTAL_READS"], row["PF_READS_ALIGNED"], row["PCT_PF_READS_ALIGNED"])
```

- **示例输出：**

```text
completed 38 38 1
```

## 7) `dnakit.sequencing_qc.gatk_collect_base_distribution_by_cycle()`

- **作用：** 绘制 SAM/BAM/CRAM 文件中每个测序周期的核苷酸分布，并输出每周期百分比。
- **处理方法：** 由带类型的参数构造 `gatk CollectBaseDistributionByCycle` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectBaseDistributionByCycle --help` 校验。 The chart output needs R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API：** `dnakit.sequencing_qc.gatk_collect_base_distribution_by_cycle(input_path: 'str', output_path: 'str', chart_output: 'str', *, aligned_reads_only: 'bool | None' = None, pf_reads_only: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `chart_output`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2852。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_base_distribution_by_cycle

work = Path(tempfile.mkdtemp())
result = gatk_collect_base_distribution_by_cycle("toy.bam", work / "base.txt", work / "base.pdf", verbose=False)
print(result["status"], (work / "base.pdf").is_file(), len(result["metrics"]) > 0)
```

- **示例输出：**

```text
completed True True
```

## 8) `dnakit.sequencing_qc.gatk_mean_quality_by_cycle()`

- **作用：** 按测序周期统计平均碱基质量并绘图。
- **处理方法：** 由带类型的参数构造 `gatk MeanQualityByCycle` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk MeanQualityByCycle --help` 校验。 The chart output needs R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API：** `dnakit.sequencing_qc.gatk_mean_quality_by_cycle(input_path: 'str', output_path: 'str', chart_output: 'str', *, aligned_reads_only: 'bool | None' = None, pf_reads_only: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `chart_output`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2883。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_mean_quality_by_cycle

work = Path(tempfile.mkdtemp())
result = gatk_mean_quality_by_cycle("toy.bam", work / "mean.txt", work / "mean.pdf", verbose=False)
print(result["status"], (work / "mean.pdf").is_file(), len(result["metrics"]) > 0)
```

- **示例输出：**

```text
completed True False
```

## 9) `dnakit.sequencing_qc.gatk_quality_score_distribution()`

- **作用：** 绘制 SAM/BAM/CRAM 文件的碱基质量值分布（每个质量值的碱基数）。
- **处理方法：** 由带类型的参数构造 `gatk QualityScoreDistribution` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk QualityScoreDistribution --help` 校验。 The chart output needs R (Rscript) with the packages of the GATK conda environment; without it GATK stops with an error.
- **API：** `dnakit.sequencing_qc.gatk_quality_score_distribution(input_path: 'str', output_path: 'str', chart_output: 'str', *, aligned_reads_only: 'bool | None' = None, pf_reads_only: 'bool | None' = None, include_no_calls: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `chart_output`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2884。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_quality_score_distribution

work = Path(tempfile.mkdtemp())
result = gatk_quality_score_distribution("toy.bam", work / "qual.txt", work / "qual.pdf", verbose=False)
print(result["status"], (work / "qual.pdf").is_file(), len(result["metrics"]) > 0)
```

- **示例输出：**

```text
completed True False
```

## 10) `dnakit.sequencing_qc.gatk_collect_hiseqx_pf_fail_metrics()`

- **作用：** 把 HiSeq X BaseCalls 目录中的 PF 失败读段分为错位、空、多克隆与未解释四类。
- **处理方法：** 由带类型的参数构造 `gatk CollectHiSeqXPfFailMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectHiSeqXPfFailMetrics --help` 校验。 Needs the BaseCalls directory of a HiSeq X run.
- **API：** `dnakit.sequencing_qc.gatk_collect_hiseqx_pf_fail_metrics(basecalls_dir: 'str', lane: 'int', output_path: 'str', *, n_cycles: 'int | None' = None, prob_explicit_reads: 'float | None' = None, num_processors: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`basecalls_dir`, `lane`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2854。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_hiseqx_pf_fail_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_hiseqx_pf_fail_metrics(basecalls_dir="missing_dir", lane=1, output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] basecalls_dir must be an existing directory
```

## 11) `dnakit.sequencing_qc.gatk_collect_hs_metrics()`

- **作用：** 收集杂交捕获（HS）指标：on-bait/on-target 碱基、富集倍数、靶区覆盖度与 GC 缺失。
- **处理方法：** 由带类型的参数构造 `gatk CollectHsMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectHsMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_hs_metrics(input_path: 'str', output_path: 'str', bait_intervals: 'Sequence[str]', target_intervals: 'Sequence[str]', *, per_target_coverage: 'str | None' = None, per_base_coverage: 'str | None' = None, theoretical_sensitivity_output: 'str | None' = None, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, bait_set_name: 'str | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_base_quality: 'int | None' = None, clip_overlapping_reads: 'bool | None' = None, coverage_cap: 'int | None' = None, near_distance: 'int | None' = None, include_indels: 'bool | None' = None, sample_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `bait_intervals`, `target_intervals`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2855。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_hs_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_hs_metrics("toy.bam", work / "hs.txt", bait_intervals=["toy.interval_list"], target_intervals=["toy.interval_list"], verbose=False)
row = result["metrics"][0]
print(row["BAIT_TERRITORY"], row["TARGET_TERRITORY"], row["ON_BAIT_BASES"])
```

- **示例输出：**

```text
350 350 860
```

## 12) `dnakit.sequencing_qc.gatk_collect_independent_replicate_metrics()`

- **作用：** 利用 VCF 中的杂合位点估计被标记为重复的读段中实为独立观测的比例（GATK 实验性工具）。
- **处理方法：** 由带类型的参数构造 `gatk CollectIndependentReplicateMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectIndependentReplicateMetrics --help` 校验。 GATK marks this tool EXPERIMENTAL (use at your own risk).
- **API：** `dnakit.sequencing_qc.gatk_collect_independent_replicate_metrics(input_path: 'str', vcf: 'str', output_path: 'str', *, matrix_output: 'str | None' = None, sample: 'str | None' = None, barcode_tag: 'str | None' = None, barcode_bq: 'str | None' = None, minimum_gq: 'int | None' = None, minimum_bq: 'int | None' = None, minimum_mq: 'int | None' = None, minimum_barcode_bq: 'int | None' = None, filter_unpaired_reads: 'bool | None' = None, reference: 'str | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `vcf`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2856。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_independent_replicate_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_independent_replicate_metrics(input_path="missing.dat", vcf="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 13) `dnakit.sequencing_qc.gatk_collect_jumping_library_metrics()`

- **作用：** 收集跳跃（mate-pair）文库指标：片段、跳跃与非跳跃配对数及其方向。
- **处理方法：** 由带类型的参数构造 `gatk CollectJumpingLibraryMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectJumpingLibraryMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_jumping_library_metrics(input_paths: 'Sequence[str]', output_path: 'str', *, chimera_kb_min: 'int | None' = None, minimum_mapping_quality: 'int | None' = None, tail_limit: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2857。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_jumping_library_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_jumping_library_metrics(["toy.bam"], work / "jump.txt", verbose=False)
row = result["metrics"][0]
print(row["JUMP_PAIRS"], row["NONJUMP_PAIRS"])
```

- **示例输出：**

```text
0 19
```

## 14) `dnakit.sequencing_qc.gatk_collect_multiple_metrics()`

- **作用：** 一次遍历文件运行多个 Picard 指标收集器（比对汇总、插入片段、周期质量、质量分布、碱基分布、GC 偏倚、测序伪迹、质量产出）。
- **处理方法：** 由带类型的参数构造 `gatk CollectMultipleMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectMultipleMetrics --help` 校验。 All collectors except CollectQualityYieldMetrics and CollectSequencingArtifactMetrics write PDF charts and need R; pass programs=['null', 'CollectQualityYieldMetrics'] to run without R.
- **API：** `dnakit.sequencing_qc.gatk_collect_multiple_metrics(input_path: 'str', output_prefix: 'str', *, reference: 'str | None' = None, programs: 'Sequence[str] | None' = None, extra_arguments: 'Sequence[str] | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, db_snp: 'str | None' = None, intervals_file: 'str | None' = None, ref_flat: 'str | None' = None, include_unpaired: 'bool | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2858。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics_files`, `collected`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_multiple_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_multiple_metrics("toy.bam", work / "multi", reference="toy.fa", programs=["null", "CollectQualityYieldMetrics"], verbose=False)
print(result["collected"])
```

- **示例输出：**

```text
['quality_yield_metrics']
```

## 15) `dnakit.sequencing_qc.gatk_collect_oxog_metrics()`

- **作用：** 收集氧化（8-氧鸟嘌呤）伪迹指标：各序列上下文中的 G>T 错误率。
- **处理方法：** 由带类型的参数构造 `gatk CollectOxoGMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectOxoGMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_oxog_metrics(input_path: 'str', output_path: 'str', reference: 'str', *, db_snp: 'str | None' = None, intervals_file: 'str | None' = None, contexts: 'Sequence[str] | None' = None, context_size: 'int | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_quality_score: 'int | None' = None, minimum_insert_size: 'int | None' = None, maximum_insert_size: 'int | None' = None, include_non_pf_reads: 'bool | None' = None, use_oq: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2859。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_oxog_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_oxog_metrics("toy.bam", work / "oxog.txt", reference="toy.fa", verbose=False)
print(len(result["metrics"]) > 0, result["metrics"][0]["SAMPLE_ALIAS"])
```

- **示例输出：**

```text
True S1
```

## 16) `dnakit.sequencing_qc.gatk_collect_quality_yield_metrics()`

- **作用：** 收集通过质量阈值与厂商过滤的读段产出指标：variant='standard'（标准）、'flow'（基于 flow 的读段文件）或 'snvq'（SNV 质量）。
- **处理方法：** 由带类型的参数构造 `gatk CollectQualityYieldMetrics / CollectQualityYieldMetricsFlow / CollectQualityYieldMetricsSNVQ` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectQualityYieldMetrics / CollectQualityYieldMetricsFlow / CollectQualityYieldMetricsSNVQ --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_quality_yield_metrics(input_path: 'str', output_path: 'str', variant: 'str', *, use_original_qualities: 'bool | None' = None, include_bq_histogram: 'bool | None' = None, alternate_quality_attribute: 'str | None' = None, include_secondary_alignments: 'bool | None' = None, include_supplemental_alignments: 'bool | None' = None, reference: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `variant`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2860, G2861, G2862。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `variant`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_quality_yield_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_quality_yield_metrics("toy.bam", work / "yield.txt", variant="standard", verbose=False)
row = result["metrics"][0]
print(result["variant"] if "variant" in result else "standard", row["TOTAL_READS"], row["TOTAL_BASES"], row["Q20_BASES"])
```

- **示例输出：**

```text
standard 38 1900 1900
```

## 17) `dnakit.sequencing_qc.gatk_collect_rnaseq_metrics()`

- **作用：** 生成 RNA-seq 比对指标：编码、UTR、内含子、基因间与核糖体区域的碱基，链特异性以及 5'/3' 偏倚。
- **处理方法：** 由带类型的参数构造 `gatk CollectRnaSeqMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectRnaSeqMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_rnaseq_metrics(input_path: 'str', output_path: 'str', ref_flat: 'str', strand_specificity: 'str', *, ribosomal_intervals: 'str | None' = None, chart_output: 'str | None' = None, end_bias_bases: 'int | None' = None, minimum_length: 'int | None' = None, rrna_fragment_percentage: 'float | None' = None, ignore_sequences: 'Sequence[str] | None' = None, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `ref_flat`, `strand_specificity`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2863。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_rnaseq_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_rnaseq_metrics("toy.bam", work / "rna.txt", ref_flat="toy.refFlat", strand_specificity="NONE", verbose=False)
row = result["metrics"][0]
print(row["PF_ALIGNED_BASES"], row["CODING_BASES"], row["INTRONIC_BASES"])
```

- **示例输出：**

```text
1900 1640 260
```

## 18) `dnakit.sequencing_qc.gatk_collect_rrbs_metrics()`

- **作用：** 收集简化代表性亚硫酸氢盐测序（RRBS）指标：CpG 转化率、覆盖度与读段起点数。
- **处理方法：** 由带类型的参数构造 `gatk CollectRrbsMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectRrbsMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_rrbs_metrics(input_path: 'str', metrics_file_prefix: 'str', reference: 'str', *, c_quality_threshold: 'int | None' = None, next_base_quality_threshold: 'int | None' = None, minimum_read_length: 'int | None' = None, max_mismatch_rate: 'float | None' = None, sequence_names: 'Sequence[str] | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, assume_sorted: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `metrics_file_prefix`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2864。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics_files`, `metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_rrbs_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_rrbs_metrics("toy.bam", work / "rrbs", reference="toy.fa", verbose=False)
print(result["status"], len(result["metrics_files"]) >= 2)
```

- **示例输出：**

```text
completed True
```

## 19) `dnakit.sequencing_qc.gatk_collect_sequencing_artifact_metrics()`

- **作用：** 量化单碱基测序伪迹（pre-adapter、bait-bias，以及按替换类型和上下文的错误汇总）。
- **处理方法：** 由带类型的参数构造 `gatk CollectSequencingArtifactMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectSequencingArtifactMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_sequencing_artifact_metrics(input_path: 'str', output_prefix: 'str', reference: 'str', *, db_snp: 'str | None' = None, intervals_file: 'str | None' = None, contexts_to_print: 'Sequence[str] | None' = None, context_size: 'int | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_quality_score: 'int | None' = None, minimum_insert_size: 'int | None' = None, maximum_insert_size: 'int | None' = None, include_duplicates: 'bool | None' = None, include_non_pf_reads: 'bool | None' = None, include_unpaired: 'bool | None' = None, tandem_reads: 'bool | None' = None, use_oq: 'bool | None' = None, file_extension: 'str | None' = None, assume_sorted: 'bool | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_prefix`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2866。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics_files`, `metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_sequencing_artifact_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_sequencing_artifact_metrics("toy.bam", work / "artifact", reference="toy.fa", verbose=False)
print(sorted(path.rsplit(".", 1)[1] for path in result["metrics_files"]))
```

- **示例输出：**

```text
['bait_bias_detail_metrics', 'bait_bias_summary_metrics', 'error_summary_metrics', 'pre_adapter_detail_metrics', 'pre_adapter_summary_metrics']
```

## 20) `dnakit.sequencing_qc.gatk_collect_targeted_pcr_metrics()`

- **作用：** 从靶向扩增子测序计算 PCR 相关指标：on-amplicon/on-target 碱基、覆盖度与灵敏度。
- **处理方法：** 由带类型的参数构造 `gatk CollectTargetedPcrMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectTargetedPcrMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_targeted_pcr_metrics(input_path: 'str', output_path: 'str', amplicon_intervals: 'str', target_intervals: 'Sequence[str]', *, per_target_coverage: 'str | None' = None, per_base_coverage: 'str | None' = None, theoretical_sensitivity_output: 'str | None' = None, custom_amplicon_set_name: 'str | None' = None, reference: 'str | None' = None, metric_accumulation_level: 'Sequence[str] | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_base_quality: 'int | None' = None, clip_overlapping_reads: 'bool | None' = None, coverage_cap: 'int | None' = None, near_distance: 'int | None' = None, include_indels: 'bool | None' = None, sample_size: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `amplicon_intervals`, `target_intervals`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2867。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_targeted_pcr_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_targeted_pcr_metrics("toy.bam", work / "pcr.txt", amplicon_intervals="toy.interval_list", target_intervals=["toy.interval_list"], verbose=False)
row = result["metrics"][0]
print(row["AMPLICON_TERRITORY"], row["ON_AMPLICON_BASES"])
```

- **示例输出：**

```text
350 860
```

## 21) `dnakit.sequencing_qc.gatk_collect_umi_prevalence_metrics()`

- **作用：** 统计重复标记 BAM 中重复集内 UMI 的个数分布（有 1、2、3 … 个不同 UMI 的重复集各有多少）。
- **处理方法：** 由带类型的参数构造 `gatk CollectUmiPrevalenceMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectUmiPrevalenceMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_umi_prevalence_metrics(input_path: 'str', output_path: 'str', *, barcode_tag: 'str | None' = None, barcode_bq: 'str | None' = None, minimum_barcode_bq: 'int | None' = None, minimum_mq: 'int | None' = None, filter_unpaired_reads: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2868。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_umi_prevalence_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_umi_prevalence_metrics(input_path="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 22) `dnakit.sequencing_qc.gatk_collect_variant_calling_metrics()`

- **作用：** 从 VCF 收集每样本与整体变异检出指标（SNP/indel 数、Ti/Tv、杂合/纯合比、dbSNP 一致性）。
- **处理方法：** 由带类型的参数构造 `gatk CollectVariantCallingMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectVariantCallingMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_variant_calling_metrics(input_path: 'str', dbsnp: 'str', output_prefix: 'str', *, target_intervals: 'str | None' = None, sequence_dictionary: 'str | None' = None, gvcf_input: 'bool | None' = None, thread_count: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `dbsnp`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2869。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics_files`, `summary_metrics`, `detail_metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_variant_calling_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_variant_calling_metrics("toy.vcf.gz", "toy.vcf.gz", work / "calling", verbose=False)
row = result["summary_metrics"][0]
print(row["TOTAL_SNPS"], row["NUM_IN_DB_SNP"], row["PCT_DBSNP"])
```

- **示例输出：**

```text
15 15 1
```

## 23) `dnakit.sequencing_qc.gatk_collect_wgs_metrics()`

- **作用：** 收集全基因组测序的覆盖度与性能指标：平均覆盖度、达到各深度阈值的碱基比例、被排除碱基比例。
- **处理方法：** 由带类型的参数构造 `gatk CollectWgsMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectWgsMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_collect_wgs_metrics(input_path: 'str', output_path: 'str', reference: 'str', *, intervals_file: 'str | None' = None, minimum_mapping_quality: 'int | None' = None, minimum_base_quality: 'int | None' = None, coverage_cap: 'int | None' = None, read_length: 'int | None' = None, count_unpaired: 'bool | None' = None, include_bq_histogram: 'bool | None' = None, use_fast_algorithm: 'bool | None' = None, sample_size: 'int | None' = None, locus_accumulation_cap: 'int | None' = None, theoretical_sensitivity_output: 'str | None' = None, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2870。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `histogram`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_collect_wgs_metrics

work = Path(tempfile.mkdtemp())
result = gatk_collect_wgs_metrics("toy.bam", work / "wgs.txt", reference="toy.fa", minimum_base_quality=1, minimum_mapping_quality=1, verbose=False)
row = result["metrics"][0]
print(row["GENOME_TERRITORY"], row["MEAN_COVERAGE"], row["PCT_EXC_DUPE"])
```

- **示例输出：**

```text
1000 1.9 0
```

## 24) `dnakit.sequencing_qc.gatk_calculate_fingerprint_metrics()`

- **作用：** 对指纹（SAM/BAM/CRAM 或 VCF）相对于单倍型图谱计算统计量，检查其可用性。
- **处理方法：** 由带类型的参数构造 `gatk CalculateFingerprintMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CalculateFingerprintMetrics --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_calculate_fingerprint_metrics(haplotype_map: 'str', input_paths: 'Sequence[str]', output_path: 'str', *, calculate_by: 'str | None' = None, genotype_lod_threshold: 'float | None' = None, number_of_sampling: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`haplotype_map`, `input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2844。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_calculate_fingerprint_metrics

work = Path(tempfile.mkdtemp())
result = gatk_calculate_fingerprint_metrics("toy.haplotype_map", ["toy.bam"], work / "fingerprint_metrics.txt", verbose=False)
row = result["metrics"][0]
print(row["SAMPLE_ALIAS"], row["HAPLOTYPES"], row["HAPLOTYPES_WITH_EVIDENCE"])
```

- **示例输出：**

```text
S1 15 13
```

## 25) `dnakit.sequencing_qc.gatk_calculate_read_group_checksum()`

- **作用：** 基于读组（RG）为 SAM/BAM/CRAM 文件生成哈希值；读组内容相同的文件哈希相同。
- **处理方法：** 由带类型的参数构造 `gatk CalculateReadGroupChecksum` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CalculateReadGroupChecksum --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_calculate_read_group_checksum(input_path: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2845。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `text`, `checksum`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_calculate_read_group_checksum

work = Path(tempfile.mkdtemp())
a = gatk_calculate_read_group_checksum("toy.bam", work / "a.txt", verbose=False)
b = gatk_calculate_read_group_checksum("umi.bam", work / "b.txt", verbose=False)
print(len(a["checksum"]), a["checksum"] == b["checksum"])
```

- **示例输出：**

```text
32 True
```

## 26) `dnakit.sequencing_qc.gatk_check_duplicate_marking()`

- **作用：** 检查同名读段（及其配对读段）的重复标记是否一致。
- **处理方法：** 由带类型的参数构造 `gatk CheckDuplicateMarking` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CheckDuplicateMarking --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_check_duplicate_marking(input_path: 'str', *, output_path: 'str | None' = None, mode: 'str | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2846。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `has_findings`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_duplicate_marking

work = Path(tempfile.mkdtemp())
result = gatk_check_duplicate_marking("toy.bam", verbose=False)
print(result["status"], result["has_findings"])
```

- **示例输出：**

```text
completed False
```

## 27) `dnakit.sequencing_qc.gatk_check_fingerprint()`

- **作用：** 从 SAM/BAM/CRAM 或 VCF 文件计算指纹，并与提供的基因型比较（同一样本假设的 LOD）。
- **处理方法：** 由带类型的参数构造 `gatk CheckFingerprint` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CheckFingerprint --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_check_fingerprint(input_path: 'str', genotypes: 'str', haplotype_map: 'str', *, output_prefix: 'str | None' = None, summary_output: 'str | None' = None, detail_output: 'str | None' = None, expected_sample_alias: 'str | None' = None, observed_sample_alias: 'str | None' = None, genotype_lod_threshold: 'float | None' = None, ignore_read_groups: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `genotypes`, `haplotype_map`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2847。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `summary_metrics`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_fingerprint

work = Path(tempfile.mkdtemp())
result = gatk_check_fingerprint("toy.vcf.gz", "toy.vcf.gz", "toy.haplotype_map", output_prefix=work / "check", verbose=False)
row = result["summary_metrics"][0]
print(row["SAMPLE"], row["HAPLOTYPES_WITH_GENOTYPES"], round(row["LOD_EXPECTED_SAMPLE"], 2))
```

- **示例输出：**

```text
S1 15 4.71
```

## 28) `dnakit.sequencing_qc.gatk_check_pileup()`

- **作用：** 把 GATK 内部的 pileup 与 samtools mpileup 的 pileup 比较并报告差异。
- **处理方法：** 由带类型的参数构造 `gatk CheckPileup` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CheckPileup --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_check_pileup(input_paths: 'Sequence[str]', pileup: 'str', reference: 'str', *, output_path: 'str | None' = None, intervals: 'Sequence[str] | None' = None, ignore_overlaps: 'bool | None' = None, continue_after_error: 'bool | None' = None, max_depth_per_sample: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `pileup`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2848。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `text`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_pileup

work = Path(tempfile.mkdtemp())
try:
    gatk_check_pileup(input_paths=["missing.dat"], pileup="missing.dat", reference="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 29) `dnakit.sequencing_qc.gatk_check_terminator_block()`

- **作用：** 检查 BGZF 压缩文件（如 BAM）的最后一个块是否为格式正确的终止块。
- **处理方法：** 由带类型的参数构造 `gatk CheckTerminatorBlock` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CheckTerminatorBlock --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_check_terminator_block(input_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2849。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`terminator_ok`, `has_findings`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_check_terminator_block

work = Path(tempfile.mkdtemp())
result = gatk_check_terminator_block("toy.bam", verbose=False)
print(result["terminator_ok"])
```

- **示例输出：**

```text
True
```

## 30) `dnakit.sequencing_qc.gatk_compare_base_qualities()`

- **作用：** 比较两个包含相同读段且顺序一致的 SAM/BAM/CRAM 文件的碱基质量。
- **处理方法：** 由带类型的参数构造 `gatk CompareBaseQualities` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CompareBaseQualities --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_compare_base_qualities(input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, throw_on_diff: 'bool | None' = None, round_down_quantized: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2871。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `text`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_compare_base_qualities

work = Path(tempfile.mkdtemp())
result = gatk_compare_base_qualities(["toy.bam", "umi.bam"], output_path=work / "bq.txt", verbose=False)
print(result["text"].splitlines()[1])
```

- **示例输出：**

```text
all 1900 quality scores are the same
```

## 31) `dnakit.sequencing_qc.gatk_compare_duplicates_spark()`

- **作用：** 判断两个可能相同的 BAM 是否具有相同的重复读段（Spark 工具，本地运行）。
- **处理方法：** 由带类型的参数构造 `gatk CompareDuplicatesSpark` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CompareDuplicatesSpark --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_compare_duplicates_spark(input_paths: 'Sequence[str]', input2: 'str', *, output_path: 'str | None' = None, output2_path: 'str | None' = None, throw_on_diff: 'bool | None' = None, print_summary: 'bool | None' = None, intervals: 'Sequence[str] | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `input2`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2872。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `summary_lines`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_compare_duplicates_spark

work = Path(tempfile.mkdtemp())
try:
    gatk_compare_duplicates_spark(input_paths=["missing.dat"], input2="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 32) `dnakit.sequencing_qc.gatk_compare_sams()`

- **作用：** 严格或宽松地比较两个 SAM/BAM/CRAM 文件（头部、比对与重复标记），并写出比较结果。
- **处理方法：** 由带类型的参数构造 `gatk CompareSAMs` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CompareSAMs --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_compare_sams(input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, lenient_dup: 'bool | None' = None, lenient_header: 'bool | None' = None, lenient_low_mq_alignment: 'bool | None' = None, lenient_unknown_mq_alignment: 'bool | None' = None, low_mq_threshold: 'int | None' = None, compare_mq: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2873。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `files_match`, `has_findings`, `text`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_compare_sams

work = Path(tempfile.mkdtemp())
result = gatk_compare_sams(["toy.bam", "umi.bam"], output_path=work / "cmp.txt", verbose=False)
print(result["files_match"], result["has_findings"])
```

- **示例输出：**

```text
True False
```

## 33) `dnakit.sequencing_qc.gatk_convert_haplotype_database_to_vcf()`

- **作用：** 把单倍型数据库（haplotype map）文件转换为 VCF。
- **处理方法：** 由带类型的参数构造 `gatk ConvertHaplotypeDatabaseToVcf` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk ConvertHaplotypeDatabaseToVcf --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_convert_haplotype_database_to_vcf(input_path: 'str', output_path: 'str', reference: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2874。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_convert_haplotype_database_to_vcf

work = Path(tempfile.mkdtemp())
result = gatk_convert_haplotype_database_to_vcf("toy.haplotype_map", work / "map.vcf", reference="toy.fa", verbose=False)
print(result["record_count"])
```

- **示例输出：**

```text
15
```

## 34) `dnakit.sequencing_qc.gatk_crosscheck_fingerprints()`

- **作用：** 检查输入文件中的所有数据是否来自同一个体（或应匹配的分组是否匹配）。
- **处理方法：** 由带类型的参数构造 `gatk CrosscheckFingerprints` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CrosscheckFingerprints --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_crosscheck_fingerprints(haplotype_map: 'str', input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, second_input_paths: 'Sequence[str] | None' = None, crosscheck_by: 'str | None' = None, crosscheck_mode: 'str | None' = None, lod_threshold: 'float | None' = None, expect_all_groups_to_match: 'bool | None' = None, output_errors_only: 'bool | None' = None, matrix_output: 'str | None' = None, num_threads: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`haplotype_map`, `input_paths`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2875。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics`, `has_findings`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_crosscheck_fingerprints

work = Path(tempfile.mkdtemp())
result = gatk_crosscheck_fingerprints("toy.haplotype_map", ["toy.bam", "umi.bam"], output_path=work / "cross.txt", verbose=False)
print(result["has_findings"], sorted({row["RESULT"] for row in result["metrics"]}))
```

- **示例输出：**

```text
False ['EXPECTED_MATCH']
```

## 35) `dnakit.sequencing_qc.gatk_dump_tabix_index()`

- **作用：** 以文本形式转储 tabix 索引文件（.tbi）：序列名、bin、chunk 与线性索引。
- **处理方法：** 由带类型的参数构造 `gatk DumpTabixIndex` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk DumpTabixIndex --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_dump_tabix_index(tabix_index: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`tabix_index`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2876。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `text`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_dump_tabix_index

work = Path(tempfile.mkdtemp())
result = gatk_dump_tabix_index("toy.vcf.gz.tbi", work / "tbi.txt", verbose=False)
print(result["text"].splitlines()[3])
```

- **示例输出：**

```text
chr1 binned index:
```

## 36) `dnakit.sequencing_qc.gatk_extract_fingerprint()`

- **作用：** 依据单倍型图谱从 SAM/BAM/CRAM 文件计算指纹并写成 VCF。
- **处理方法：** 由带类型的参数构造 `gatk ExtractFingerprint` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk ExtractFingerprint --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_extract_fingerprint(input_path: 'str', haplotype_map: 'str', output_path: 'str', reference: 'str', *, contamination: 'float | None' = None, sample_alias: 'str | None' = None, locus_max_reads: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `haplotype_map`, `output_path`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2878。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_extract_fingerprint

work = Path(tempfile.mkdtemp())
result = gatk_extract_fingerprint("toy.bam", "toy.haplotype_map", work / "fingerprint.vcf", reference="toy.fa", verbose=False)
print(result["record_count"])
```

- **示例输出：**

```text
15
```

## 37) `dnakit.sequencing_qc.gatk_get_sample_name()`

- **作用：** 把 SAM/BAM/CRAM 文件中唯一的样本名写入文件。
- **处理方法：** 由带类型的参数构造 `gatk GetSampleName` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GetSampleName --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_get_sample_name(input_paths: 'Sequence[str]', output_path: 'str', *, use_url_encoding: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2880。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `sample_name`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_get_sample_name

work = Path(tempfile.mkdtemp())
result = gatk_get_sample_name(["toy.bam"], work / "sample.txt", verbose=False)
print(result["sample_name"])
```

- **示例输出：**

```text
S1
```

## 38) `dnakit.sequencing_qc.gatk_identify_contaminant()`

- **作用：** 在给定污染估计的情况下，从 SAM/BAM/CRAM 文件计算污染样本的指纹并写成 VCF。
- **处理方法：** 由带类型的参数构造 `gatk IdentifyContaminant` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk IdentifyContaminant --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_identify_contaminant(input_path: 'str', haplotype_map: 'str', output_path: 'str', reference: 'str', *, contamination: 'float | None' = None, extract_contaminated: 'bool | None' = None, sample_alias: 'str | None' = None, locus_max_reads: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `haplotype_map`, `output_path`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2881。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_identify_contaminant

work = Path(tempfile.mkdtemp())
result = gatk_identify_contaminant("toy.bam", "toy.haplotype_map", work / "contaminant.vcf", reference="toy.fa", contamination=0.1, verbose=False)
print(result["record_count"])
```

- **示例输出：**

```text
15
```

## 39) `dnakit.sequencing_qc.gatk_lift_over_haplotype_map()`

- **作用：** 用 chain 文件把单倍型数据库（haplotype map）从一个参考基因组版本转换到另一个版本。
- **处理方法：** 由带类型的参数构造 `gatk LiftOverHaplotypeMap` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk LiftOverHaplotypeMap --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_lift_over_haplotype_map(input_path: 'str', chain: 'str', sequence_dictionary: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `chain`, `sequence_dictionary`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2882。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_lift_over_haplotype_map

work = Path(tempfile.mkdtemp())
result = gatk_lift_over_haplotype_map("old.haplotype_map", "shift.chain", "toy.dict", work / "lifted.haplotype_map", verbose=False)
print((work / "lifted.haplotype_map").read_text().splitlines()[-1].split("\t")[:3])
```

- **示例输出：**

```text
['chr1', '40', 'rs1']
```

## 40) `dnakit.sequencing_qc.gatk_validate_sam_file()`

- **作用：** 校验 SAM/BAM/CRAM 文件（头部、标志、配对、CIGAR、索引），并报告错误与警告。
- **处理方法：** 由带类型的参数构造 `gatk ValidateSamFile` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk ValidateSamFile --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_validate_sam_file(input_path: 'str', *, mode: 'str | None' = None, output_path: 'str | None' = None, ignore: 'Sequence[str] | None' = None, ignore_warnings: 'bool | None' = None, max_output: 'int | None' = None, is_bisulfite_sequenced: 'bool | None' = None, skip_mate_validation: 'bool | None' = None, index_validation_stringency: 'str | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2885。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `valid`, `has_findings`, `problems`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_validate_sam_file

work = Path(tempfile.mkdtemp())
result = gatk_validate_sam_file("toy.bam", mode="SUMMARY", verbose=False)
print(result["valid"], result["has_findings"])
```

- **示例输出：**

```text
True False
```

## 41) `dnakit.sequencing_qc.gatk_view_sam()`

- **作用：** 输出 SAM/BAM/CRAM 文件的记录（及头部），可只选已比对/未比对、PF/非 PF 读段。
- **处理方法：** 由带类型的参数构造 `gatk ViewSam` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk ViewSam --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_view_sam(input_path: 'str', alignment_status: 'str', pf_status: 'str', *, header_only: 'bool | None' = None, records_only: 'bool | None' = None, interval_list: 'str | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `alignment_status`, `pf_status`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2886。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`header_lines`, `record_count`, `records`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_view_sam

work = Path(tempfile.mkdtemp())
result = gatk_view_sam("toy.bam", alignment_status="All", pf_status="All", records_only=True, verbose=False)
print(result["record_count"], result["records"][0].split("\t")[:6])
```

- **示例输出：**

```text
38 ['r1', '99', 'chr1', '1', '60', '50M']
```

## 42) `dnakit.sequencing_qc.gatk_htsget_download()`

- **作用：** 通过 GA4GH htsget 端点（服务器 URL 上的记录 ID）下载文件到本地；需要网络。
- **处理方法：** 由带类型的参数构造 `gatk HtsgetReader` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk HtsgetReader --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_htsget_download(url: 'str', record_id: 'str', output_path: 'str', *, data_format: 'str | None' = None, data_class: 'str | None' = None, intervals: 'str | None' = None, check_md5: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`url`, `record_id`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2887。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_htsget_download

work = Path(tempfile.mkdtemp())
try:
    gatk_htsget_download(url="ftp://x", record_id="x", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] url must be an http(s) URL of an htsget endpoint.
```

## 43) `dnakit.sequencing_qc.gatk_split_cram()`

- **作用：** 按指定记录数把 CRAM 文件高效拆分为较小的 CRAM 文件（不重新压缩容器）。
- **处理方法：** 由带类型的参数构造 `gatk SplitCRAM` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk SplitCRAM --help` 校验。
- **API：** `dnakit.sequencing_qc.gatk_split_cram(input_path: 'str', output_dir: 'str', *, shard_records: 'int | None' = None, shard_max_output_count: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_dir`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2895。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `shards`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequencing_qc.gatk_p00_metrics_cli_functions import gatk_split_cram

work = Path(tempfile.mkdtemp())
result = gatk_split_cram("toy.cram", work / "shards", shard_records=10, verbose=False)
print([path.rsplit("/", 1)[1] for path in result["shards"]])
```

- **示例输出：**

```text
['shard_0000.cram', 'shard_0001.cram']
```
