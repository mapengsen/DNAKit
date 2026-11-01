# GATK 命令行封装：区间列表、基因分型芯片、flow 读段与 PathSeq

对 `gatk` 命令行的带类型封装，涵盖 Picard 区间列表（转换、比较、集合运算、坐标转换、分散）、Illumina 基因分型芯片（GTC 转 VCF、清单、ADPC、芯片指标）、Ultima 式 flow 读段工具以及 PathSeq 参考构建。 每个函数构造固定的 `gatk` 命令行（不经 shell），只写入新路径，并返回可 JSON 序列化的字典；GATK 采用 BSD-3-Clause 许可（Picard 为 MIT），不随包分发。函数位于 `dnakit.annotation.gatk_p00_intervals_cli_functions`，并从 `dnakit.annotation` 导出；示例从模块导入，因此两种情况下都能运行。

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

## 1) `dnakit.annotation.gatk_convert_interval_list()`

- **作用：** 在 BED 与 Picard interval list 格式之间转换：output_format='interval_list' 对应 BedToIntervalList，output_format='bed' 对应 IntervalListToBed。
- **处理方法：** 由带类型的参数构造 `gatk BedToIntervalList / IntervalListToBed` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk BedToIntervalList / IntervalListToBed --help` 校验。
- **API：** `dnakit.annotation.gatk_convert_interval_list(input_path: 'str', output_path: 'str', output_format: 'str', *, sequence_dictionary: 'str | None' = None, keep_length_zero_intervals: 'bool | None' = None, unique: 'bool | None' = None, score: 'int | None' = None, sort: 'bool | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `output_format`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2903, G2905。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `output_format`, `interval_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_convert_interval_list

work = Path(tempfile.mkdtemp())
to_list = gatk_convert_interval_list("toy.bed", work / "from_bed.interval_list", output_format="interval_list", sequence_dictionary="toy.dict", verbose=False)
to_bed = gatk_convert_interval_list("toy.interval_list", work / "from_list.bed", output_format="bed", verbose=False)
print(to_list["interval_count"], (work / "from_list.bed").read_text().splitlines()[0])
```

- **示例输出：**

```text
2 chr1	100	300	t1	500	+
```

## 2) `dnakit.annotation.gatk_compare_interval_lists()`

- **作用：** 比较两个 interval list 是否相等，并报告不同的区间。
- **处理方法：** 由带类型的参数构造 `gatk CompareIntervalLists` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CompareIntervalLists --help` 校验。
- **API：** `dnakit.annotation.gatk_compare_interval_lists(first: 'str', second: 'str', reference: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`first`, `second`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2904。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`identical`, `has_findings`, `report_lines`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_compare_interval_lists

work = Path(tempfile.mkdtemp())
same = gatk_compare_interval_lists("toy.interval_list", "toy.interval_list", reference="toy.fa", verbose=False)
print(same["identical"], same["report_lines"])
```

- **示例输出：**

```text
True ['Intervals are equal']
```

## 3) `dnakit.annotation.gatk_interval_list_tools()`

- **作用：** 操作 Picard interval list：拼接、并集、交集、差集、对称差、填充、反转、合并、按倍数切分以及分散。
- **处理方法：** 由带类型的参数构造 `gatk IntervalListTools` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk IntervalListTools --help` 校验。
- **API：** `dnakit.annotation.gatk_interval_list_tools(input_paths: 'Sequence[str]', *, output_path: 'str | None' = None, output_directory: 'str | None' = None, action: 'str | None' = None, second_input_paths: 'Sequence[str] | None' = None, padding: 'int | None' = None, invert: 'bool | None' = None, unique: 'bool | None' = None, sort: 'bool | None' = None, dont_merge_abutting: 'bool | None' = None, break_bands_at_multiples_of: 'int | None' = None, scatter_count: 'int | None' = None, scatter_content: 'int | None' = None, subdivision_mode: 'str | None' = None, count_output: 'str | None' = None, output_value: 'str | None' = None, comments: 'Sequence[str] | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2906。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `interval_count`, `scattered_lists`, `count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_interval_list_tools

work = Path(tempfile.mkdtemp())
padded = gatk_interval_list_tools(["toy.interval_list"], output_path=work / "padded.interval_list", padding=10, verbose=False)
print((work / "padded.interval_list").read_text().splitlines()[-2:])
```

- **示例输出：**

```text
['chr1\t91\t310\t+\tt1', 'chr2\t41\t210\t+\tt2']
```

## 4) `dnakit.annotation.gatk_lift_over_interval_list()`

- **作用：** 用 chain 文件把 interval list 从一个参考基因组版本转换到另一个版本；无法转换的区间被拒绝。
- **处理方法：** 由带类型的参数构造 `gatk LiftOverIntervalList` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk LiftOverIntervalList --help` 校验。
- **API：** `dnakit.annotation.gatk_lift_over_interval_list(input_path: 'str', chain: 'str', sequence_dictionary: 'str', output_path: 'str', *, reject: 'str | None' = None, min_liftover_pct: 'float | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `chain`, `sequence_dictionary`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2907。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `interval_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_lift_over_interval_list

work = Path(tempfile.mkdtemp())
result = gatk_lift_over_interval_list("old.interval_list", "shift.chain", "toy.dict", work / "lifted.interval_list", reject=work / "rejected.interval_list", verbose=False)
print((work / "lifted.interval_list").read_text().splitlines()[-1])
```

- **示例输出：**

```text
chr1	111	210	+	t1
```

## 5) `dnakit.annotation.gatk_split_intervals()`

- **作用：** 把 interval list（或整个参考序列）拆成多个子区间文件，用于分散-汇聚并行处理。
- **处理方法：** 由带类型的参数构造 `gatk SplitIntervals` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk SplitIntervals --help` 校验。
- **API：** `dnakit.annotation.gatk_split_intervals(output_dir: 'str', *, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, scatter_count: 'int | None' = None, subdivision_mode: 'str | None' = None, dont_mix_contigs: 'bool | None' = None, min_contig_size: 'int | None' = None, interval_file_prefix: 'str | None' = None, extension: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`output_dir`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2909。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `interval_files`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_split_intervals

work = Path(tempfile.mkdtemp())
result = gatk_split_intervals(work / "scatter", intervals=["toy.interval_list"], reference="toy.fa", scatter_count=2, verbose=False)
print([path.rsplit("/", 1)[1] for path in result["interval_files"]])
```

- **示例输出：**

```text
['0000-scattered.interval_list', '0001-scattered.interval_list']
```

## 6) `dnakit.annotation.gatk_bpm_to_normalization_manifest_csv()`

- **作用：** 把 Illumina 微珠池清单（.bpm）与聚类文件转换为 bpm.csv 归一化清单。
- **处理方法：** 由带类型的参数构造 `gatk BpmToNormalizationManifestCsv` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk BpmToNormalizationManifestCsv --help` 校验。
- **API：** `dnakit.annotation.gatk_bpm_to_normalization_manifest_csv(input_path: 'str', cluster_file: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `cluster_file`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2896。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `row_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_bpm_to_normalization_manifest_csv

work = Path(tempfile.mkdtemp())
try:
    gatk_bpm_to_normalization_manifest_csv(input_path="missing.dat", cluster_file="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 7) `dnakit.annotation.gatk_combine_genotyping_array_vcfs()`

- **作用：** 把多个基因分型芯片 VCF 合并为一个 VCF。
- **处理方法：** 由带类型的参数构造 `gatk CombineGenotypingArrayVcfs` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CombineGenotypingArrayVcfs --help` 校验。
- **API：** `dnakit.annotation.gatk_combine_genotyping_array_vcfs(input_paths: 'Sequence[str]', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2897。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_combine_genotyping_array_vcfs

work = Path(tempfile.mkdtemp())
try:
    gatk_combine_genotyping_array_vcfs(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 8) `dnakit.annotation.gatk_compare_gtc_files()`

- **作用：** 对照微珠池清单比较两个 Illumina GTC 文件（基因型与信号强度）。
- **处理方法：** 由带类型的参数构造 `gatk CompareGtcFiles` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CompareGtcFiles --help` 校验。
- **API：** `dnakit.annotation.gatk_compare_gtc_files(input_paths: 'Sequence[str]', illumina_bead_pool_manifest_file: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `illumina_bead_pool_manifest_file`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2898。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`files_match`, `has_findings`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_compare_gtc_files

work = Path(tempfile.mkdtemp())
try:
    gatk_compare_gtc_files(input_paths=["missing.dat"], illumina_bead_pool_manifest_file="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must hold exactly two GTC files.
```

## 9) `dnakit.annotation.gatk_create_extended_illumina_manifest()`

- **作用：** 创建 GtcToVcf 所需的扩展 Illumina 清单（含参考基因组特异的映射位置与等位基因）。
- **处理方法：** 由带类型的参数构造 `gatk CreateExtendedIlluminaManifest` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CreateExtendedIlluminaManifest --help` 校验。
- **API：** `dnakit.annotation.gatk_create_extended_illumina_manifest(input_path: 'str', output_path: 'str', reference_sequence: 'str', report_file: 'str', *, bad_assays_file: 'str | None' = None, cluster_file: 'str | None' = None, dbsnp_file: 'str | None' = None, target_build: 'str | None' = None, supported_builds: 'Sequence[str] | None' = None, supported_reference_files: 'Sequence[str] | None' = None, supported_chain_files: 'Sequence[str] | None' = None, flag_duplicates: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `reference_sequence`, `report_file`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2899。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `row_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_create_extended_illumina_manifest

work = Path(tempfile.mkdtemp())
try:
    gatk_create_extended_illumina_manifest(input_path="missing.dat", output_path=work / "o", reference_sequence="missing.dat", report_file=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 10) `dnakit.annotation.gatk_gtc_to_vcf()`

- **作用：** 借助扩展清单、微珠池清单与聚类文件，把 Illumina GTC 基因型文件转换为 VCF。
- **处理方法：** 由带类型的参数构造 `gatk GtcToVcf` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GtcToVcf --help` 校验。
- **API：** `dnakit.annotation.gatk_gtc_to_vcf(input_path: 'str', output_path: 'str', extended_illumina_manifest: 'str', illumina_bead_pool_manifest_file: 'str', cluster_file: 'str', sample_alias: 'str', reference: 'str', *, expected_gender: 'str | None' = None, gender_gtc: 'str | None' = None, fingerprint_genotypes_vcf_file: 'str | None' = None, analysis_version_number: 'int | None' = None, pipeline_version: 'str | None' = None, do_not_allow_calls_on_zeroed_out_assays: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `output_path`, `extended_illumina_manifest`, `illumina_bead_pool_manifest_file`, `cluster_file`, `sample_alias`, `reference`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2900。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_gtc_to_vcf

work = Path(tempfile.mkdtemp())
try:
    gatk_gtc_to_vcf(input_path="missing.dat", output_path=work / "o", extended_illumina_manifest="missing.dat", illumina_bead_pool_manifest_file="missing.dat", cluster_file="missing.dat", sample_alias="x", reference="missing.dat", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 11) `dnakit.annotation.gatk_merge_ped_into_vcf()`

- **作用：** 把 zCall 产生的单样本 PED 文件合并进单样本 VCF。
- **处理方法：** 由带类型的参数构造 `gatk MergePedIntoVcf` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk MergePedIntoVcf --help` 校验。
- **API：** `dnakit.annotation.gatk_merge_ped_into_vcf(original_vcf: 'str', ped_file: 'str', map_file: 'str', zcall_thresholds_file: 'str', zcall_version: 'str', output_path: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`original_vcf`, `ped_file`, `map_file`, `zcall_thresholds_file`, `zcall_version`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2901。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_merge_ped_into_vcf

work = Path(tempfile.mkdtemp())
try:
    gatk_merge_ped_into_vcf(original_vcf="missing.dat", ped_file="missing.dat", map_file="missing.dat", zcall_thresholds_file="missing.dat", zcall_version="x", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] original_vcf must be an existing file
```

## 12) `dnakit.annotation.gatk_vcf_to_adpc()`

- **作用：** 把基因分型芯片 VCF 转换为用于 verifyIDintensity 的 ADPC（adpc.bin）文件，并输出样本名与标记数文件。
- **处理方法：** 由带类型的参数构造 `gatk VcfToAdpc` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk VcfToAdpc --help` 校验。
- **API：** `dnakit.annotation.gatk_vcf_to_adpc(vcf_paths: 'Sequence[str]', output_path: 'str', samples_file: 'str', num_markers_file: 'str', *, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`vcf_paths`, `output_path`, `samples_file`, `num_markers_file`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2902。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `num_markers`, `samples`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_vcf_to_adpc

work = Path(tempfile.mkdtemp())
try:
    gatk_vcf_to_adpc(vcf_paths=["missing.dat"], output_path=work / "o", samples_file=work / "o", num_markers_file=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] vcf_paths must be an existing file
```

## 13) `dnakit.annotation.gatk_collect_arrays_variant_calling_metrics()`

- **作用：** 从基因分型芯片 VCF 收集每样本（detail）与整体（summary）指标：检出率、纯合/杂合数、dbSNP 一致性。
- **处理方法：** 由带类型的参数构造 `gatk CollectArraysVariantCallingMetrics` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CollectArraysVariantCallingMetrics --help` 校验。
- **API：** `dnakit.annotation.gatk_collect_arrays_variant_calling_metrics(input_path: 'str', dbsnp: 'str', output_prefix: 'str', *, call_rate_pf_threshold: 'float | None' = None, sequence_dictionary: 'str | None' = None, num_processors: 'int | None' = None, reference: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_path`, `dbsnp`, `output_prefix`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2851。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `metrics_files`, `summary_metrics`, `detail_metrics`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_collect_arrays_variant_calling_metrics

work = Path(tempfile.mkdtemp())
try:
    gatk_collect_arrays_variant_calling_metrics(input_path="missing.dat", dbsnp="missing.dat", output_prefix=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_path must be an existing file
```

## 14) `dnakit.annotation.gatk_add_flow_base_quality()`

- **作用：** 给基于 flow 的 SAM/BAM/CRAM 文件的读段添加碱基质量属性（Ultima 式 flow 数据）。
- **处理方法：** 由带类型的参数构造 `gatk AddFlowBaseQuality` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk AddFlowBaseQuality --help` 校验。 Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API：** `dnakit.annotation.gatk_add_flow_base_quality(input_paths: 'Sequence[str]', output_path: 'str', *, maximal_quality_score: 'int | None' = None, minimal_error_rate: 'float | None' = None, replace_quality_mode: 'bool | None' = None, keep_boundary_flows: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2888。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_add_flow_base_quality

work = Path(tempfile.mkdtemp())
try:
    gatk_add_flow_base_quality(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 15) `dnakit.annotation.gatk_add_flow_snv_quality()`

- **作用：** 给基于 flow 的 CRAM/BAM 的读段添加 SNV 质量属性（Ultima 式 flow 数据）。
- **处理方法：** 由带类型的参数构造 `gatk AddFlowSNVQuality` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk AddFlowSNVQuality --help` 校验。 Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API：** `dnakit.annotation.gatk_add_flow_snv_quality(input_paths: 'Sequence[str]', output_path: 'str', *, snvq_mode: 'str | None' = None, output_quality_attribute: 'str | None' = None, max_phred_score: 'float | None' = None, keep_supplementary_alignments: 'bool | None' = None, include_qc_failed_read: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2889。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_add_flow_snv_quality

work = Path(tempfile.mkdtemp())
try:
    gatk_add_flow_snv_quality(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 16) `dnakit.annotation.gatk_calculate_average_combined_annotations()`

- **作用：** 把 GenomicsDB 对样本求和的注释除以样本数，得到平均值。
- **处理方法：** 由带类型的参数构造 `gatk CalculateAverageCombinedAnnotations` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk CalculateAverageCombinedAnnotations --help` 校验。
- **API：** `dnakit.annotation.gatk_calculate_average_combined_annotations(variant: 'str', summed_annotations_to_divide: 'Sequence[str]', output_path: 'str', *, intervals: 'Sequence[str] | None' = None, reference: 'str | None' = None, variant_output_filtering: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`variant`, `summed_annotations_to_divide`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2890。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_calculate_average_combined_annotations

work = Path(tempfile.mkdtemp())
try:
    gatk_calculate_average_combined_annotations(variant="missing.dat", summed_annotations_to_divide=["x"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] variant must be an existing file
```

## 17) `dnakit.annotation.gatk_flow_feature_mapper()`

- **作用：** 在基于 flow 的 BAM 中映射特征（目前为 SNV），并写成 VCF。
- **处理方法：** 由带类型的参数构造 `gatk FlowFeatureMapper` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk FlowFeatureMapper --help` 校验。 Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API：** `dnakit.annotation.gatk_flow_feature_mapper(input_paths: 'Sequence[str]', output_path: 'str', *, mapping_feature: 'str | None' = None, copy_attrs: 'Sequence[str] | None' = None, copy_attr_prefix: 'str | None' = None, min_score: 'float | None' = None, max_score: 'float | None' = None, limit_score: 'float | None' = None, snv_identical_bases: 'int | None' = None, snv_identical_bases_after: 'int | None' = None, include_dup_reads: 'bool | None' = None, include_qc_failed_reads: 'bool | None' = None, keep_negatives: 'bool | None' = None, exclude_nan_scores: 'bool | None' = None, keep_supplementary_alignments: 'bool | None' = None, emit_ref_confidence: 'str | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2891。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `record_count`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_flow_feature_mapper

work = Path(tempfile.mkdtemp())
try:
    gatk_flow_feature_mapper(input_paths=["missing.dat"], output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 18) `dnakit.annotation.gatk_flow_pair_hmm_align_reads_to_haplotypes()`

- **作用：** 用基于 flow 的 PairHMM 生成 读段 x 单倍型 矩阵，给出每条读段对每个单倍型的对数似然。
- **处理方法：** 由带类型的参数构造 `gatk FlowPairHMMAlignReadsToHaplotypes` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk FlowPairHMMAlignReadsToHaplotypes --help` 校验。 Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API：** `dnakit.annotation.gatk_flow_pair_hmm_align_reads_to_haplotypes(input_paths: 'Sequence[str]', haplotypes: 'str', output_path: 'str', *, aligner: 'str | None' = None, likelihood_calculation_engine: 'str | None' = None, ref_haplotype: 'str | None' = None, concise_output_format: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `haplotypes`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2892。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `matrix_written`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_flow_pair_hmm_align_reads_to_haplotypes

work = Path(tempfile.mkdtemp())
try:
    gatk_flow_pair_hmm_align_reads_to_haplotypes(input_paths=["missing.dat"], haplotypes="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 19) `dnakit.annotation.gatk_ground_truth_reads_builder()`

- **作用：** 利用读段与母本、父本单倍型参考序列，生成用于碱基识别训练的灵活稳健的真值读段集（CSV）。
- **处理方法：** 由带类型的参数构造 `gatk GroundTruthReadsBuilder` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GroundTruthReadsBuilder --help` 校验。 Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API：** `dnakit.annotation.gatk_ground_truth_reads_builder(input_paths: 'Sequence[str]', maternal_ref: 'str', paternal_ref: 'str', ancestral_translators_base_path: 'str', output_csv: 'str', *, subsampling_ratio: 'float | None' = None, max_output_reads: 'int | None' = None, min_haplotype_score: 'float | None' = None, min_mq: 'float | None' = None, include_supp_align: 'bool | None' = None, fill_trimmed_reads: 'bool | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `maternal_ref`, `paternal_ref`, `ancestral_translators_base_path`, `output_csv`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2893。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_ground_truth_reads_builder

work = Path(tempfile.mkdtemp())
try:
    gatk_ground_truth_reads_builder(input_paths=["missing.dat"], maternal_ref="missing.dat", paternal_ref="missing.dat", ancestral_translators_base_path=work / "o", output_csv=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 20) `dnakit.annotation.gatk_ground_truth_scorer()`

- **作用：** 把基于 flow 的读段与参考序列或真值比较评分，输出逐读段 CSV 和可选报告。
- **处理方法：** 由带类型的参数构造 `gatk GroundTruthScorer` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk GroundTruthScorer --help` 校验。 Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK marks the flow tools experimental.
- **API：** `dnakit.annotation.gatk_ground_truth_scorer(input_paths: 'Sequence[str]', output_csv: 'str', *, report_file: 'str | None' = None, features_file: 'str | None' = None, genome_prior: 'str | None' = None, add_mean_call: 'bool | None' = None, exclude_zero_flows: 'bool | None' = None, omit_zeros_from_report: 'bool | None' = None, quality_percentiles: 'str | None' = None, use_softclipped_bases: 'bool | None' = None, normalized_score_threshold: 'float | None' = None, reference: 'str | None' = None, intervals: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`input_paths`, `output_csv`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2894。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_ground_truth_scorer

work = Path(tempfile.mkdtemp())
try:
    gatk_ground_truth_scorer(input_paths=["missing.dat"], output_csv=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] input_paths must be an existing file
```

## 21) `dnakit.annotation.gatk_pathseq_build_kmers()`

- **作用：** 构建宿主参考序列的 k-mer 集合（哈希集合；给出假阳性概率时为 Bloom 过滤器），供 PathSeq 扣除宿主读段使用。
- **处理方法：** 由带类型的参数构造 `gatk PathSeqBuildKmers` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk PathSeqBuildKmers --help` 校验。
- **API：** `dnakit.annotation.gatk_pathseq_build_kmers(reference: 'str', output_path: 'str', *, kmer_size: 'int | None' = None, kmer_mask: 'str | None' = None, kmer_spacing: 'int | None' = None, bloom_false_positive_probability: 'float | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`reference`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2910。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`, `kmer_files`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_pathseq_build_kmers

work = Path(tempfile.mkdtemp())
result = gatk_pathseq_build_kmers("toy.fa", work / "host", kmer_size=15, verbose=False)
print([path.rsplit("/", 1)[1] for path in result["kmer_files"]])
```

- **示例输出：**

```text
['host.hss']
```

## 22) `dnakit.annotation.gatk_pathseq_build_reference_taxonomy()`

- **作用：** 用 NCBI 分类学数据转储与参考 FASTA 构建 PathSeq 所需的微生物参考分类数据文件（.db）。
- **处理方法：** 由带类型的参数构造 `gatk PathSeqBuildReferenceTaxonomy` 的固定参数列表（不经 shell），在输出文件所在目录运行并记录版本，要求退出码为 0 且预期输出文件存在。该工具的其他长选项通过 `options` 传入，并按 `gatk PathSeqBuildReferenceTaxonomy --help` 校验。
- **API：** `dnakit.annotation.gatk_pathseq_build_reference_taxonomy(reference: 'str', tax_dump: 'str', output_path: 'str', *, refseq_catalog: 'str | None' = None, genbank_catalog: 'str | None' = None, min_non_virus_contig_length: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **依赖：** 需要 PATH 中的 `gatk` 启动脚本（GATK 4.6.2.0 或更高，Java 17），或设置 `DNAKIT_GATK` / 传入 `gatk_path`；缺失时抛出带安装提示的 `BackendUnavailableError`。
- **输入：** 必填参数：`reference`, `tax_dump`, `output_path`。所有可选参数默认为 None，即沿用 GATK 默认值。输出路径必须尚不存在。对应工作清单行 G2911。
- **输出：** 标准结果字典（`status`、`tool`、`version`、`command`、`return_code`、`outputs`、`log_tail`、`provenance`），另含：`outputs`。
- **示例说明：** 示例数据没有适用于该工具的有效输入（它需要无法缩小的数据，如 Illumina 运行目录、GTC 文件、flow 读段，或 Python gcnvkernel 环境），因此示例展示的是启动 GATK 之前执行的输入检查。

- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.annotation.gatk_p00_intervals_cli_functions import gatk_pathseq_build_reference_taxonomy

work = Path(tempfile.mkdtemp())
try:
    gatk_pathseq_build_reference_taxonomy(reference="missing.dat", tax_dump="missing.dat", output_path=work / "o", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] reference must be an existing file
```
