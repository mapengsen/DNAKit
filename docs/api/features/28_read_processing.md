# 测序数据整理

## 1) 匹配双端 reads {#bio-10}

- **作用：** 按名称匹配 R1/R2，并识别缺少配对的记录。
- **计算方法：** 去掉 ID 末尾的 `/1`、`/2` 后匹配，缺失的一端返回 `None`。
- **API：** `dnakit.bioinformatics.reads.pair_reads(first, second, max_pending=100000, include_orphans=True)`
- **输入：** 两个 `Read` 迭代器；每端规范化后的 ID 不能重复。
- **输出：** 返回 `ReadPair` 迭代器，每项包含配对 ID 和两端 reads。
- **示例代码：**

```python
from dnakit.bioinformatics.reads import Read, pair_reads

pairs = pair_reads([Read("r1/1", "ACGT")], [Read("r1/2", "TGCA")])
print([(pair.id, pair.first.sequence, pair.second.sequence) for pair in pairs])
```

- **示例结果：**

```text
[('r1', 'ACGT', 'TGCA')]
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 6) 长读长过滤与统计 {#bio-long-read-qc}

DNAKit 为 NanoFilt、Filtlong 和 NanoPlot 提供统一的有界外部接口。调用者
需要单独安装对应程序；DNAKit 不复制 GPL/MIT 后端，也不下载数据库。

- **NanoFilt API：** `dnakit.sequencing_qc.filter_nanopore_reads(source, output_directory, min_quality=None, min_length=None, headcrop=None, tailcrop=None, ...)`。
- **Filtlong API：** `dnakit.sequencing_qc.filter_long_reads(source, output_directory, min_length=None, keep_percent=None, target_bases=None, ...)`。
- **NanoPlot API：** `dnakit.sequencing_qc.summarize_long_reads(source, output_directory, threads=1, min_length=None, min_quality=None, ...)`。
- **Porechop API：** `dnakit.sequencing_qc.trim_porechop(source, output_directory, threads=1, no_split=False, ...)`。
- **Trimmomatic API：** `dnakit.sequencing_qc.trim_trimmomatic(read1, output_directory, read2=None, adapter_fasta=None, leading=None, trailing=None, min_length=None, ...)`。
- **AfterQC API：** `dnakit.sequencing_qc.run_afterqc(read1, output_directory, read2=None, qc_only=False, ...)`。
- **输入与输出：** 两个过滤接口读取本地 FASTQ 并捕获 `filtered.fastq` 和 stderr；NanoPlot 使用 FASTQ 模式、关闭静态图片，并检查 `NanoStats.txt` 与 `NanoPlot-report.html`。输出目录必须是新的目录，所有路径、参数、超时和产物大小都会校验。
- **Agent/MCP：** 三个函数登记为 `external_process`，调用时仍需显式允许外部副作用；接口只覆盖固定的质量、长度、裁剪、保留比例、目标碱基数和线程参数。
- **引用与许可：** NanoFilt 引用 De Coster et al. (2018), [10.1093/bioinformatics/bty149](https://doi.org/10.1093/bioinformatics/bty149)，GPL-3.0；Filtlong 官方项目未提供 DOI，GPL-3.0；NanoPlot 引用 De Coster & Rademakers (2023), [10.1093/bioinformatics/btad311](https://doi.org/10.1093/bioinformatics/btad311)，MIT。NanoFilt 官方仓库已说明项目不再更新。
- **Porechop 引用与许可：** Wick et al. (2017), [10.1099/mgen.0.000132](https://doi.org/10.1099/mgen.0.000132)，GPL-3.0；官方仓库已说明项目不再维护。
- **Trimmomatic 引用与许可：** Bolger, Lohse & Usadel (2014), [10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170)，GPL-3.0；adapter 序列由 Illumina 另行提供，使用时要遵守上游条款。
- **AfterQC 引用与许可：** Chen et al. (2017), [10.1186/s12859-017-1469-3](https://doi.org/10.1186/s12859-017-1469-3)，MIT。AfterQC 官方 README 建议维护中的 C++ 工具 fastp；此接口保留 AfterQC 论文工作流的可复现入口。

这些门面不会宣称覆盖长读长工具的全部命令行；需要 Filtlong 外部参考、单位后缀或 NanoPlot 其他输入格式时，应直接按上游文档运行对应程序，并在论文中保留上游引用。

## 2) 恢复损坏 FASTQ 中的合法记录 {#bio-11}

- **作用：** 从损坏文件中保留可识别的四行 FASTQ 记录。
- **计算方法：** 扫描记录头、序列、分隔符与等长质量行，并记录丢弃行数。
- **API：** `dnakit.bioinformatics.reads.repair_fastq(source, destination, phred_offset=33)`
- **输入：** 损坏 FASTQ 和新输出路径；多行 FASTQ 应用 `read_fastx` 严格读取。
- **输出：** 返回 `RepairResult`，包含修复后 FASTQ 路径、恢复的 reads 数和丢弃的行数。
- **示例代码：**

```python
from tempfile import TemporaryDirectory
from pathlib import Path
from dnakit.bioinformatics.reads import repair_fastq

with TemporaryDirectory() as folder:
    source = Path(folder) / "broken.fq"
    source.write_text("invalid\n@r1\nACGT\n+\nIIII\n")
    result = repair_fastq(source, Path(folder) / "repaired.fq")
    print(result.recovered_reads, result.discarded_lines)
```

- **示例结果：**

```text
1 1
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 3) 排序、取交集、重命名和校验 reads {#bio-12}

- **作用：** 整理记录顺序、提取共有序列并计算序列摘要。
- **计算方法：** 排序或固定种子打乱；按 ID/序列取交集；正则重命名；按记录边界计算 SHA256。
- **API：** `dnakit.bioinformatics.reads.arrange_reads(reads, operation='sort')`; `dnakit.bioinformatics.reads.common_reads(collections)`; `dnakit.bioinformatics.reads.rename_reads(reads, pattern, replacement)`; `dnakit.bioinformatics.reads.sequence_digest(reads)`
- **输入：** `Read` 集合；摘要与记录顺序有关，与 ID 无关。
- **输出：** 返回排序、交集或重命名后的 `Read` 记录；校验函数返回比较结果。
- **示例代码：**

```python
from dnakit.bioinformatics.reads import Read, arrange_reads, common_reads, rename_reads, sequence_digest

reads = [Read("b", "TGCA"), Read("a", "ACGT")]
ordered = arrange_reads(reads)
shared = common_reads([ordered, [Read("a", "ACGT")]])
renamed = list(rename_reads(shared, "^a$", "sample_a"))
print([read.id for read in ordered])
print([read.id for read in renamed])
print(sequence_digest(shared) == sequence_digest(renamed))
```

- **示例结果：**

```text
['a', 'b']
['sample_a']
True
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 4) 整理、标记重复并索引 BAM {#bio-16}

- **作用：** 把 SAM/BAM 整理为可用于变异检测的 BAM。
- **计算方法：** samtools 依次执行 name-sort、fixmate、坐标排序、markdup、index 和 stats。
- **API：** `dnakit.bioinformatics.ngs.prepare_bam(source, output_directory, mark_duplicates=True, config=None)`
- **输入：** 已有 `aligned.sam` 或 BAM；samtools 在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
- **Agent/MCP API：** `dnakit.ngs.prepare_bam(source, output_directory, mark_duplicates=True, ...)`；返回 JSON 可序列化的产物路径、后端和日志。
- **引用与许可：** SAMtools 使用 Danecek et al. (2021), DOI `10.1093/gigascience/giab008`（MIT）；后端由调用者安装。
- **示例代码：**

```python
from dnakit.bioinformatics.ngs import prepare_bam

result = prepare_bam("aligned.sam", "new_bam")
print(result.artifacts["aligned.bam"])
```

- **示例结果：**

返回 `ToolResult`，含 `aligned.bam`、索引和统计文件。外部程序示例未执行。

[依赖与实现说明](../../faq/bioinformatics.md#reads)。

## 5) 外部测序 QC 工具 {#bio-external-qc}

DNAKit 为 Cutadapt、fastp、FastQC 和 MultiQC 提供同一套 Python/Agent/MCP 外部后端边界。它们是实际可执行的有界适配器：调用者必须自己安装后端，DNAKit 不捆绑可执行程序、Java、插件或数据库。

- **共同方法：** DNAKit 只生成固定参数白名单的 argv，使用无 shell 的有界执行器；会检查可执行文件、输入路径、超时、输出大小、退出码和声明的产物，并拒绝覆盖已有输出。任意未列入白名单的命令行参数都会被拒绝。
- **共同输出：** 成功时返回 JSON 可序列化信封，包含 `status`、`backend`、实际 argv、`outputs`、截断后的后端输出和 `provenance`。缺少后端、非零退出、缺失/空报告、无效 JSON 报告和资源限制错误都会以结构化异常返回。
- **Agent/MCP：** 四个函数在工具目录中标记为 `external_process`；MCP 调用仍需调用方显式允许外部副作用。该接口覆盖 DNAKit 的固定参数子集，不声称重导出上游完整 CLI。

### Cutadapt

- **API：** `dnakit.sequencing_qc.run_cutadapt(input, *, options=None, verbose=True)`。
- **输入：** `input` 是 FASTQ 路径，或 `{"read1": path, "read2": path}`。`options` 需要 `executable_path`、`output`/`output1`；双端输入还需要 `output2`，并至少提供一个接头参数，例如 `adapters=["AGATCGGAAG"]`。支持固定的 `cores`、`minimum_overlap`、`minimum_length`、`error_rate`、`quality_cutoff`、`discard_untrimmed` 和 `json_report` 等选项。
- **输出：** 一个或两个 trimmed FASTQ，以及可选 JSON 报告；所有输出必须位于已有目录且不能覆盖已有文件。
- **引用与许可：** Martin (2011), [10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)；Cutadapt 上游 MIT 许可。

### fastp

- **API：** `dnakit.sequencing_qc.run_fastp(input, *, options=None, verbose=True)`。
- **输入：** 单端或 `{"read1": path, "read2": path}` 双端 FASTQ。`options` 需要 `executable_path`、`output`/`output1`，双端还需要 `output2`；可选质量、长度、线程、接头和过滤开关。
- **输出：** 一个或两个过滤后的 FASTQ，以及 JSON 和 HTML 报告。若未指定报告路径，DNAKit 会在输出 FASTQ 同目录生成不冲突的默认报告名。
- **引用与许可：** Chen (2025), [10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078)；fastp 上游 MIT 许可。

### FastQC

- **API：** `dnakit.sequencing_qc.run_fastqc(input, *, options=None, verbose=True)`。
- **输入：** 一个到八个本地 FASTQ/BAM/SAM 风格输入路径；`options` 需要 `executable_path` 和已有的空 `output_directory`。适配器按 FastQC 0.12.x 接口支持 `threads`、`noextract`、`nogroup`、`quiet`、`format`、`contaminant`、`adapter` 和 2–10 的 `kmer`。
- **输出：** 每个输入对应的 `<sample>_fastqc.html` 与 `<sample>_fastqc.zip`，并检查两个文件确实由后端生成。一次调用最多监控八个输入（十六个报告），以保持外部产物限制。
- **引用与许可：** [FastQC 官方项目页](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)；调查表未提供 DOI；上游 GPL-3.0 许可。当前参数契约针对 FastQC 0.12.x，其他版本应先单独核对 CLI。

### MultiQC

- **API：** `dnakit.sequencing_qc.run_multiqc(input, *, options=None, verbose=True)`。
- **输入：** 一个或多个本地文件/目录；`options` 需要 `executable_path` 和已有的空 `output_directory`，并以 `filename` 指定 `.html` 报告名。支持 `title`、`interactive`、`zip_data`、`flat`、`dirs`、`dirs_depth` 和本地 `config` 等固定选项；`zip_data` 映射到 MultiQC 的 `--zip-data-dir`。
- **输出：** 指定的 HTML 汇总报告。空输出目录让 DNAKit 能明确检查报告是否生成，避免把后端“退出码为 0 但没有结果”误报为成功。
- **引用与许可：** Ewels et al. (2016), [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)；MultiQC 上游 GPL-3.0-or-later 许可。

**示例代码（Cutadapt）：**

```python
from pathlib import Path
from dnakit.sequencing_qc import run_cutadapt

result = run_cutadapt(
    {"read1": Path("reads_R1.fastq"), "read2": Path("reads_R2.fastq")},
    options={
        "executable_path": "/usr/local/bin/cutadapt",
        "output1": Path("results/trimmed_R1.fastq"),
        "output2": Path("results/trimmed_R2.fastq"),
        "adapters": ["AGATCGGAAG"],
        "json_report": Path("results/cutadapt.json"),
    },
    verbose=False,
)
print(result["outputs"])
```

[依赖与实现说明](../../faq/bioinformatics.md#reads)。
