# GATK 比对与读段文件工具（GATK，P01 包）

这组 35 个函数位于 `dnakit.alignment.gatk_p01_alignment_cli_functions`，调用 Broad Institute 的 [GATK](https://github.com/broadinstitute/gatk) 4.7.0.0 命令行（含其内置的 Picard 工具）。每个函数运行一个固定工具，或由模式参数选择的一族工具；不经过 shell，带超时和输出大小限制。`options` 只接受 `gatk <Tool> --help` 列出的长选项，已由参数设置的选项不能重复，GATK 会执行脚本的 `JAVASCRIPT_FILE` 不可用。输出写到显式的新路径，不会覆盖已有文件。

GATK 是外部程序：请安装（例如 `conda install -c bioconda gatk4` 或发布包）并准备 Java 17，然后把 `gatk` 启动器放进 PATH、设置 `DNAKIT_GATK` 或传入 `gatk_path`。找不到时抛出带安装提示的 `BackendUnavailableError`（代码 `GATK_UNAVAILABLE`）；工具运行失败抛出 `BackendExecutionError`（代码 `GATK_EXECUTION_FAILED`，附日志尾部）；工具返回 0 却未写出所请求的输出时抛出 `GATK_OUTPUT_MISSING`。DNAKit 不捆绑 GATK。Spark 工具以本地模式运行（可通过 `options` 设置 `spark-master`）。参考 FASTA 旁需要 `.fai` 和 `.dict`，BAM/VCF 需要各自索引，与 GATK 本身要求一致。结果为字典，每个函数的小节列出其额外字段。

## 1) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_pathseq()`

- **作用:** 运行 PathSeq 病原发现流程：宿主读段过滤、病原参考比对、丰度评分，或一次完成全部步骤。
- **处理方法:** 按 `step` 运行 `gatk PathSeqFilterSpark`、`PathSeqBwaSpark`、`PathSeqScoreSpark` 或 `PathSeqPipelineSpark`（Spark 本地模式）。宿主和病原资源（k-mer 文件、宿主 BWA 镜像等）通过 `filter-bwa-image`、`kmer-file` 等 GATK 选项传入。 `step`: `filter` = 宿主读段过滤，`align` = 病原比对，`score` = 丰度评分，`pipeline` = 全部步骤。
- **API:** `gatk_pathseq(*, step: 'str' = 'filter', input_path: 'str | None' = None, paired_input: 'str | None' = None, unpaired_input: 'str | None' = None, microbe_bwa_image: 'str | None' = None, microbe_dict: 'str | None' = None, taxonomy_file: 'str | None' = None, paired_output: 'str | None' = None, unpaired_output: 'str | None' = None, scores_output: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：PathSeqBwaSpark, PathSeqFilterSpark, PathSeqPipelineSpark, PathSeqScoreSpark。
- **输入:** `input_path` (部分模式必需): 待检测的 BAM/CRAM/SAM 读段（宿主过滤的输入）。; `paired_input` (可选): 过滤步骤产生的配对读段 BAM。; `unpaired_input` (可选): 过滤步骤产生的未配对读段 BAM。; `microbe_bwa_image` (部分模式必需): 病原参考序列的 BWA-MEM 索引镜像（可用 `gatk_bwa_mem_index_image` 生成）。; `microbe_dict` (部分模式必需): 病原参考序列的序列字典。; `taxonomy_file` (部分模式必需): 由 PathSeqBuildReferenceTaxonomy 构建的分类数据库。; `paired_output` (可选): 配对读段输出 BAM。; `unpaired_output` (可选): 未配对读段输出 BAM。; `scores_output` (部分模式必需): 病原丰度评分表（制表符分隔）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段、`mode`，以及所请求的输出路径。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_pathseq

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_pathseq(
            input_path=str(input_path),
            paired_output=str(folder / "paired_output.out"),
            step='filter',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 2) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_detect_cram_issue_8768()`

- **作用:** 检查 CRAM 文件是否存在 GATK issue 8768 造成的碱基损坏。
- **处理方法:** 运行 `gatk CRAMIssue8768Detector`；错配率阈值和 TSV 输出通过 GATK 选项（`mismatch-rate-threshold`、`output-tsv`）设置。
- **API:** `gatk_detect_cram_issue_8768(cram_path: 'str', reference_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CRAMIssue8768Detector。
- **输入:** `cram_path` (必需): 待检查的 CRAM 文件。; `reference_path` (必需): 该 CRAM 使用的参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的报告文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为报告。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_detect_cram_issue_8768

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    cram_path = folder / "cram_path.dat"; cram_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_detect_cram_issue_8768(
            cram_path=str(cram_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 3) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_create_bam_splitting_index()`

- **作用:** 为 BAM 文件创建 Hadoop-BAM 切分索引（.sbi）。
- **处理方法:** 运行 `gatk CreateHadoopBamSplittingIndex`；索引 `<bam>.sbi` 写在 BAM 旁边，且事先不能存在。
- **API:** `gatk_create_bam_splitting_index(bam_path: 'str', *, create_bai: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CreateHadoopBamSplittingIndex。
- **输入:** `bam_path` (必需): 按坐标排序的 BAM 文件。; `create_bai` (可选): 同时写出 .bai 索引。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.sbi` 为切分索引路径。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_create_bam_splitting_index

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    bam_path = folder / "bam_path.dat"; bam_path.write_bytes(b"x")
    try:
        gatk_create_bam_splitting_index(
            bam_path=str(bam_path),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 4) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_print_bgzf_block_information()`

- **作用:** 报告 BGZF 文件每个压缩块的偏移、压缩大小和解压大小。
- **处理方法:** 运行 `gatk PrintBGZFBlockInformation` 并解析其 `Block #n at file offset ...` 列表（返回前 10000 个块，`block_count` 为总数）。
- **API:** `gatk_print_bgzf_block_information(bgzf_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：PrintBGZFBlockInformation。
- **输入:** `bgzf_path` (必需): BGZF 压缩文件（BAM、bgzip 的 VCF 等）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段，加上 `block_count` 和 `blocks`。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_print_bgzf_block_information

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    bgzf_path = folder / "bgzf_path.dat"; bgzf_path.write_bytes(b"x")
    try:
        gatk_print_bgzf_block_information(
            bgzf_path=str(bgzf_path),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 5) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_print_file_diagnostics()`

- **作用:** 写出基因组学文件的诊断信息（记录数、头信息和格式细节）。
- **处理方法:** 运行 `gatk PrintFileDiagnostics`；`count-limit` 限制检查的记录数。
- **API:** `gatk_print_file_diagnostics(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：PrintFileDiagnostics。
- **输入:** `input_path` (必需): 工具可分析格式的基因组学文件（CRAM；GATK 4.7.0.0 不支持 BAM）。; `output_path` (必需): 新建的本地诊断文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为诊断文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_print_file_diagnostics

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_print_file_diagnostics(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 6) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_read_anonymizer()`

- **作用:** 把读段的碱基替换为其比对位置的参考碱基，从而去除个人序列信息。
- **处理方法:** 运行 `gatk ReadAnonymizer`。
- **API:** `gatk_read_anonymizer(input_path: 'str', reference_path: 'str', output_path: 'str', *, ref_base_quality: 'int | None' = None, use_simple_cigar: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ReadAnonymizer。
- **输入:** `input_path` (必需): 输入 BAM/CRAM/SAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的输出比对文件。; `ref_base_quality` (可选): 替换用参考碱基的质量值。; `use_simple_cigar` (可选): 把 CIGAR 改写为整条读段的简单匹配。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为匿名化后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_read_anonymizer

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_read_anonymizer(
            input_path=str(input_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 7) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_add_comments_to_bam()`

- **作用:** 给 SAM/BAM/CRAM 文件的头添加注释（@CO）行。
- **处理方法:** 运行基于 Picard 的 `gatk AddCommentsToBam`；`comments` 的每一项对应一个 `--COMMENT`。
- **API:** `gatk_add_comments_to_bam(input_path: 'str', output_path: 'str', comments: 'Sequence[str]', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：AddCommentsToBam。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建的输出文件。; `comments` (必需): 注释文本；每条写成一行 @CO 头。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为新文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_add_comments_to_bam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_add_comments_to_bam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            comments=["processed by dnakit"],
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 8) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_add_original_alignment_tags()`

- **作用:** 把每条读段的当前比对信息记入 OA 标签（以及原始配偶 contig）。
- **处理方法:** 运行 Picard `AddOATag`（`mode=oa_tag`）或 GATK `AddOriginalAlignmentTags`（`mode=original_alignment`，同时添加原始配偶 contig 标签）。 `mode`: `oa_tag` = AddOATag，`original_alignment` = AddOriginalAlignmentTags。
- **API:** `gatk_add_original_alignment_tags(input_path: 'str', output_path: 'str', *, mode: 'str' = 'oa_tag', interval_list: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：AddOATag, AddOriginalAlignmentTags。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建的输出文件。; `interval_list` (可选): 只有与这些区间重叠的读段才会被加 OA 标签。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`outputs.output_path` 为新文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_add_original_alignment_tags

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_add_original_alignment_tags(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            mode='oa_tag',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 9) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_base_recalibrator()`

- **作用:** 生成碱基质量值重校准（BQSR）表。
- **处理方法:** 运行 `gatk BaseRecalibrator`；协变量和质量设置通过 `low-quality-tail`、`maximum-cycle-value` 等 GATK 选项给出。
- **API:** `gatk_base_recalibrator(input_path: 'str', reference_path: 'str', known_sites: 'Sequence[str]', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：BaseRecalibrator。
- **输入:** `input_path` (必需): 输入 BAM/CRAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `known_sites` (必需): 用于屏蔽的已知变异 VCF 文件（需已索引）。; `output_path` (必需): 新建的重校准表。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为重校准表。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_base_recalibrator

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    known_sites = folder / "known_sites.dat"; known_sites.write_bytes(b"x")
    try:
        gatk_base_recalibrator(
            input_path=str(input_path),
            reference_path=str(reference_path),
            known_sites=[str(known_sites)],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 10) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_apply_bqsr()`

- **作用:** 把碱基质量值重校准表应用到读段上。
- **处理方法:** 运行 `gatk ApplyBQSR`。
- **API:** `gatk_apply_bqsr(input_path: 'str', recalibration_table: 'str', output_path: 'str', *, reference_path: 'str | None' = None, static_quantized_quals: 'Sequence[str] | None' = None, use_original_qualities: 'bool | None' = None, emit_original_quals: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ApplyBQSR。
- **输入:** `input_path` (必需): 输入 BAM/CRAM。; `recalibration_table` (必需): `gatk_base_recalibrator` 生成的重校准表。; `output_path` (必需): 新建的重校准 BAM/CRAM。; `reference_path` (可选): 参考 FASTA（CRAM 需要）。; `static_quantized_quals` (可选): 量化到的固定质量值水平。; `use_original_qualities` (可选): 使用 OQ 标签作为输入质量值。; `emit_original_quals` (可选): 把原始质量值保存在 OQ 标签。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为重校准后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_apply_bqsr

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    recalibration_table = folder / "recalibration_table.dat"; recalibration_table.write_bytes(b"x")
    try:
        gatk_apply_bqsr(
            input_path=str(input_path),
            recalibration_table=str(recalibration_table),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 11) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_bqsr_pipeline_spark()`

- **作用:** 在一次 Spark 并行运行中完成 BQSR 的两个步骤（BaseRecalibrator 与 ApplyBQSR）。
- **处理方法:** 运行 `gatk BQSRPipelineSpark`（beta）；Spark 主节点（默认 `local[*]`）通过 `options` 传入 `spark-master`。
- **API:** `gatk_bqsr_pipeline_spark(input_path: 'str', reference_path: 'str', known_sites: 'Sequence[str]', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：BQSRPipelineSpark。
- **输入:** `input_path` (必需): 输入 BAM/CRAM。; `reference_path` (必需): 参考 FASTA（Spark 也支持 .2bit）。; `known_sites` (必需): 已知变异 VCF 文件（需已索引）。; `output_path` (必需): 新建的重校准 BAM。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为重校准后的 BAM。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_bqsr_pipeline_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    known_sites = folder / "known_sites.dat"; known_sites.write_bytes(b"x")
    try:
        gatk_bqsr_pipeline_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            known_sites=[str(known_sites)],
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 12) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_bam_to_bfq()`

- **作用:** 把 BAM 文件转换为 BFQ（二进制 FASTQ）文件，用于 Maq 比对器。
- **处理方法:** 运行 Picard `gatk BamToBfq`；写入 `analysis_dir` 的文件列在 `files` 中。
- **API:** `gatk_bam_to_bfq(input_path: 'str', analysis_dir: 'str', paired_run: 'bool', *, output_file_prefix: 'str | None' = None, flowcell_barcode: 'str | None' = None, lane: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：BamToBfq。
- **输入:** `input_path` (必需): 输入 BAM/SAM。; `analysis_dir` (必需): 存放 .bfq 文件的新目录。; `output_file_prefix` (可选): .bfq 文件名前缀（与 `flowcell_barcode` 二选一）。; `flowcell_barcode` (可选): 用于生成文件名的流动槽条形码（与 `output_file_prefix` 二选一）。; `paired_run` (必需): 是否为双端测序。; `lane` (可选): 要写出的泳道号。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `files`（`analysis_dir` 内的文件名）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_bam_to_bfq

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_bam_to_bfq(
            input_path=str(input_path),
            analysis_dir=str(folder / "analysis_dir"),
            output_file_prefix="reads",
            paired_run=False,
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 13) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_bwa_spark()`

- **作用:** 在 Spark 上用 BWA-MEM 把读段比对到参考序列，可选择随后标记重复。
- **处理方法:** 运行 `gatk BwaSpark`（`pipeline=align`）或 `BwaAndMarkDuplicatesPipelineSpark`（`pipeline=align_mark_duplicates`，输入需按名称排序）；两者均为 beta 工具。 `pipeline`: `align` = BwaSpark，`align_mark_duplicates` = BwaAndMarkDuplicatesPipelineSpark。
- **API:** `gatk_bwa_spark(input_path: 'str', reference_path: 'str', output_path: 'str', *, pipeline: 'str' = 'align', bwa_mem_index_image: 'str | None' = None, single_end_alignment: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：BwaAndMarkDuplicatesPipelineSpark, BwaSpark。
- **输入:** `input_path` (必需): 未比对（或按名称排序）的 BAM/CRAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的比对 BAM。; `bwa_mem_index_image` (可选): 预先生成的 BWA-MEM 索引镜像（见 `gatk_bwa_mem_index_image`）。; `single_end_alignment` (可选): 按单端比对。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`outputs.output_path` 为比对 BAM。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_bwa_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_bwa_spark(
            input_path=str(input_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            pipeline='align',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 14) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_clean_sam()`

- **作用:** 清理 SAM/BAM/CRAM 文件：对伸出参考序列末端的比对做软剪切，并把未比对读段的 MAPQ 设为 0。
- **处理方法:** 运行 Picard `gatk CleanSam`。
- **API:** `gatk_clean_sam(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CleanSam。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建的清理后文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为清理后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_clean_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_clean_sam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 15) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_clip_reads()`

- **作用:** 按序列、循环或质量剪切读段，并报告剪切统计。
- **处理方法:** 运行 `gatk ClipReads`；`clip_representation` 取 GATK 的 `WRITE_NS`（默认）、`WRITE_Q0S`、`WRITE_NS_Q0S`、`SOFTCLIP_BASES`、`HARDCLIP_BASES` 或 `REVERT_SOFTCLIPPED_BASES`。
- **API:** `gatk_clip_reads(input_path: 'str', output_path: 'str', *, clip_sequences: 'Sequence[str] | None' = None, cycles_to_trim: 'str | None' = None, q_trimming_threshold: 'int | None' = None, clip_representation: 'str | None' = None, statistics_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ClipReads。
- **输入:** `input_path` (必需): 输入 BAM/CRAM/SAM。; `output_path` (必需): 新建的剪切后文件。; `clip_sequences` (可选): 要在读段中剪掉的序列。; `cycles_to_trim` (可选): 要修剪的循环区间，例如 `1-5,10-12`。; `q_trimming_threshold` (可选): 修剪质量低于此值的末端。; `clip_representation` (可选): 被剪切碱基的表示方式。; `statistics_path` (可选): 新建的剪切统计文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为剪切后的文件（若请求则含 `statistics_path`）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_clip_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_clip_reads(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 16) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_convert_headerless_bam_shard()`

- **作用:** 给无头 BAM 分片（Hadoop-BAM 输出）加上头，转换为可读的 BAM。
- **处理方法:** 运行 `gatk ConvertHeaderlessHadoopBamShardToBam`（beta）。
- **API:** `gatk_convert_headerless_bam_shard(shard_path: 'str', header_bam_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ConvertHeaderlessHadoopBamShardToBam。
- **输入:** `shard_path` (必需): Spark 工具写出的无头 BAM 分片。; `header_bam_path` (必需): 提供头信息的 BAM 文件。; `output_path` (必需): 新建的可读 BAM。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 BAM。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_convert_headerless_bam_shard

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    shard_path = folder / "shard_path.dat"; shard_path.write_bytes(b"x")
    header_bam_path = folder / "header_bam_path.dat"; header_bam_path.write_bytes(b"x")
    try:
        gatk_convert_headerless_bam_shard(
            shard_path=str(shard_path),
            header_bam_path=str(header_bam_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 17) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_downsample_by_duplicate_set()`

- **作用:** 按给定比例丢弃 UMI 分组 BAM 中的重复集合。
- **处理方法:** 运行 `gatk DownsampleByDuplicateSet`（beta）。
- **API:** `gatk_downsample_by_duplicate_set(input_path: 'str', output_path: 'str', *, fraction_to_keep: 'float | None' = None, min_reads: 'int | None' = None, min_per_strand_reads: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：DownsampleByDuplicateSet。
- **输入:** `input_path` (必需): 带 MI 标签（`<编号>/A`、`<编号>/B`）并按 MI 排序的 UMI 分组 BAM。; `output_path` (必需): 新建的降采样 BAM。; `fraction_to_keep` (可选): 保留的重复集合比例。; `min_reads` (可选): 集合至少包含的读段数。; `min_per_strand_reads` (可选): 集合中每条链至少的读段数。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为降采样后的 BAM。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_downsample_by_duplicate_set

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_downsample_by_duplicate_set(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 18) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_downsample_sam()`

- **作用:** 随机保留 SAM/BAM/CRAM 文件中一定比例的读段。
- **处理方法:** 运行 Picard `gatk DownsampleSam`；策略默认取 GATK 自身默认值（`ConstantMemory`）。
- **API:** `gatk_downsample_sam(input_path: 'str', output_path: 'str', *, probability: 'float | None' = None, strategy: 'str | None' = None, accuracy: 'float | None' = None, random_seed: 'int | None' = None, metrics_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：DownsampleSam。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建的降采样文件。; `probability` (可选): 保留每对读段的概率（0-1）。; `strategy` (可选): 降采样策略。; `accuracy` (可选): 保留比例的精度（ConstantMemory 不适用）。; `random_seed` (可选): 随机种子（默认 1）。; `metrics_path` (可选): 新建的指标文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为降采样后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_downsample_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_downsample_sam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 19) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_position_based_downsample_sam()`

- **作用:** 按读段在流动槽芯片上的物理位置降采样（保留整片相邻簇，适用于重复研究）。
- **处理方法:** 运行 Picard `gatk PositionBasedDownsampleSam`。
- **API:** `gatk_position_based_downsample_sam(input_path: 'str', output_path: 'str', fraction: 'float', *, stop_after: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：PositionBasedDownsampleSam。
- **输入:** `input_path` (必需): 按坐标排序、读段名为 Illumina 格式的输入。; `output_path` (必需): 新建的降采样文件。; `fraction` (必需): 保留的读段位置（芯片区域）比例（0-1）。; `stop_after` (可选): 读取这么多条读段后停止（0 为全部）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为降采样后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_position_based_downsample_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_position_based_downsample_sam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            fraction=0.5,
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 20) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_extract_original_alignment_records_by_name()`

- **作用:** 按文件中列出的读段名取 BAM 子集（保留其原始比对）。
- **处理方法:** 运行 `gatk ExtractOriginalAlignmentRecordsByNameSpark`（beta）。
- **API:** `gatk_extract_original_alignment_records_by_name(input_path: 'str', read_name_file: 'str', output_path: 'str', *, invert_match: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ExtractOriginalAlignmentRecordsByNameSpark。
- **输入:** `input_path` (必需): 输入 BAM。; `read_name_file` (必需): 每行一个读段名的文本文件。; `output_path` (必需): 新建的、含所选读段的 BAM。; `invert_match` (可选): 保留未列出的读段。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为子集。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_extract_original_alignment_records_by_name

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    read_name_file = folder / "read_name_file.dat"; read_name_file.write_bytes(b"x")
    try:
        gatk_extract_original_alignment_records_by_name(
            input_path=str(input_path),
            read_name_file=str(read_name_file),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 21) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_filter_sam_reads()`

- **作用:** 按比对状态、读段名、区间或标签值取 SAM/BAM/CRAM 文件的读段子集。
- **处理方法:** 运行 Picard `gatk FilterSamReads`。不提供 `includeJavascript` 过滤，因为它会执行用户脚本。
- **API:** `gatk_filter_sam_reads(input_path: 'str', output_path: 'str', filter: 'str', *, read_list_file: 'str | None' = None, interval_list: 'str | None' = None, tag: 'str | None' = None, tag_values: 'Sequence[str] | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：FilterSamReads。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建的过滤后文件。; `filter` (必需): 使用的过滤方式。; `read_list_file` (可选): 读段名列表（read-list 过滤使用）。; `interval_list` (可选): `includePairedIntervals` 使用的区间。; `tag` (可选): 标签值过滤使用的标签名。; `tag_values` (可选): 标签值过滤使用的标签值。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为过滤后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_filter_sam_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_filter_sam_reads(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            filter="includeAligned",
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 22) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_fix_misencoded_base_quality_reads()`

- **作用:** 修正以旧版 Phred+64 编码存储的 Illumina 碱基质量值。
- **处理方法:** 运行 `gatk FixMisencodedBaseQualityReads`。
- **API:** `gatk_fix_misencoded_base_quality_reads(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：FixMisencodedBaseQualityReads。
- **输入:** `input_path` (必需): 质量值为 Illumina 编码（Phred+64）的输入。; `output_path` (必需): 新建的、质量值为 Phred+33 的文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为修正后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_fix_misencoded_base_quality_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_fix_misencoded_base_quality_reads(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 23) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_left_align_indels()`

- **作用:** 对 SAM/BAM/CRAM 文件读段 CIGAR 中的插入缺失做左对齐。
- **处理方法:** 运行 `gatk LeftAlignIndels`。
- **API:** `gatk_left_align_indels(input_path: 'str', reference_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：LeftAlignIndels。
- **输入:** `input_path` (必需): 输入 BAM/CRAM/SAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的、插入缺失已左对齐的文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为新文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_left_align_indels

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_left_align_indels(
            input_path=str(input_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 24) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_mark_duplicates_with_mate_cigar()`

- **作用:** 利用配偶 CIGAR 标记重复读段（标准、简化或 UMI 感知）。
- **处理方法:** 按 `variant` 运行 Picard 的 `MarkDuplicatesWithMateCigar`、`SimpleMarkDuplicatesWithMateCigar` 或 `UmiAwareMarkDuplicatesWithMateCigar`；指标文件解析到 `tables.metrics_path`。 `variant`: `standard` = MarkDuplicatesWithMateCigar，`simple` = SimpleMarkDuplicatesWithMateCigar（实验性），`umi_aware` = UmiAwareMarkDuplicatesWithMateCigar（实验性）。
- **API:** `gatk_mark_duplicates_with_mate_cigar(input_path: 'str', output_path: 'str', metrics_path: 'str', *, variant: 'str' = 'standard', umi_metrics_path: 'str | None' = None, remove_duplicates: 'bool | None' = None, umi_tag_name: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：MarkDuplicatesWithMateCigar, SimpleMarkDuplicatesWithMateCigar, UmiAwareMarkDuplicatesWithMateCigar。
- **输入:** `input_path` (必需): 按坐标排序、带 MC（配偶 CIGAR）标签的 BAM。; `output_path` (必需): 新建的、已标记重复的 BAM。; `metrics_path` (必需): 新建的重复指标文件。; `umi_metrics_path` (部分模式必需): 新建的 UMI 指标文件。; `remove_duplicates` (可选): 直接移除重复而不是仅标记。; `umi_tag_name` (可选): 保存 UMI 的读段标签（默认 RX）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段、`mode` 和 `tables.metrics_path`（Picard 指标行）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_mark_duplicates_with_mate_cigar

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_mark_duplicates_with_mate_cigar(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            metrics_path=str(folder / "metrics_path.out"),
            variant='standard',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 25) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_merge_bam_alignment()`

- **作用:** 把比对器输出 BAM 的比对信息与未比对 BAM 的原始数据合并。
- **处理方法:** 运行 Picard `gatk MergeBamAlignment`。
- **API:** `gatk_merge_bam_alignment(unmapped_bam: 'str', reference_path: 'str', output_path: 'str', *, aligned_bams: 'Sequence[str] | None' = None, paired_run: 'bool | None' = None, sort_order: 'str | None' = None, add_mate_cigar: 'bool | None' = None, clip_adapters: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：MergeBamAlignment。
- **输入:** `unmapped_bam` (必需): 原始未比对 BAM（按查询名排序）。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的合并 BAM。; `aligned_bams` (可选): 比对器输出（按查询名排序）；省略则只输出未比对读段。; `paired_run` (可选): 是否为双端测序。; `sort_order` (可选): 输出的排序方式。; `add_mate_cigar` (可选): 添加 MC 标签。; `clip_adapters` (可选): 剪切 MarkIlluminaAdapters 标记的接头。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为合并后的 BAM。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_merge_bam_alignment

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    unmapped_bam = folder / "unmapped_bam.dat"; unmapped_bam.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_merge_bam_alignment(
            unmapped_bam=str(unmapped_bam),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 26) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_print_distant_mates()`

- **作用:** 把配偶比对在远处的读段设为未比对并写出。
- **处理方法:** 运行 `gatk PrintDistantMates`。
- **API:** `gatk_print_distant_mates(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：PrintDistantMates。
- **输入:** `input_path` (必需): 输入 BAM/CRAM/SAM。; `output_path` (必需): 新建的文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为新文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_print_distant_mates

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_print_distant_mates(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 27) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_reorder_sam()`

- **作用:** 按另一个参考序列的 contig 顺序重排 SAM/BAM/CRAM 文件的读段。
- **处理方法:** 运行 Picard `gatk ReorderSam`。
- **API:** `gatk_reorder_sam(input_path: 'str', output_path: 'str', sequence_dictionary: 'str', *, allow_incomplete_dict_concordance: 'bool | None' = None, allow_contig_length_discordance: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ReorderSam。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建的重排后文件。; `sequence_dictionary` (必需): 目标 contig 顺序所取自的 SAM/BAM/CRAM/dict/fasta。; `allow_incomplete_dict_concordance` (可选): 允许目标中缺少部分 contig。; `allow_contig_length_discordance` (可选): 允许 contig 长度不一致。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为重排后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_reorder_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    sequence_dictionary = folder / "sequence_dictionary.dat"; sequence_dictionary.write_bytes(b"x")
    try:
        gatk_reorder_sam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            sequence_dictionary=str(sequence_dictionary),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 28) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_revert_base_quality_scores()`

- **作用:** 把碱基质量值还原为 OQ 标签中保存的原始值（可同时添加配偶 CIGAR 标签）。
- **处理方法:** 运行 `gatk RevertBaseQualityScores`（`mode=original_qualities`）或 Picard `RevertOriginalBaseQualitiesAndAddMateCigar`（`mode=with_mate_cigar`）。 `mode`: `original_qualities` = RevertBaseQualityScores，`with_mate_cigar` = RevertOriginalBaseQualitiesAndAddMateCigar。
- **API:** `gatk_revert_base_quality_scores(input_path: 'str', output_path: 'str', *, mode: 'str' = 'original_qualities', restore_original_qualities: 'bool | None' = None, sort_order: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：RevertBaseQualityScores, RevertOriginalBaseQualitiesAndAddMateCigar。
- **输入:** `input_path` (必需): 带 OQ 标签的输入 BAM/CRAM/SAM。; `output_path` (必需): 新建的文件。; `restore_original_qualities` (可选): 还原 OQ 值。; `sort_order` (可选): 输出排序方式。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`outputs.output_path` 为新文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_revert_base_quality_scores

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_revert_base_quality_scores(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            mode='original_qualities',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 29) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_revert_sam()`

- **作用:** 把 SAM/BAM/CRAM 文件还原为较早状态（未比对、原始质量值、无重复标记），以便重新比对。
- **处理方法:** 运行 Picard `gatk RevertSam`，输出为单个文件；本函数不提供 `OUTPUT_BY_READGROUP`。
- **API:** `gatk_revert_sam(input_path: 'str', output_path: 'str', *, restore_original_qualities: 'bool | None' = None, remove_duplicate_information: 'bool | None' = None, remove_alignment_information: 'bool | None' = None, sanitize: 'bool | None' = None, sample_alias: 'str | None' = None, library_name: 'str | None' = None, sort_order: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：RevertSam。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建的还原后文件。; `restore_original_qualities` (可选): 还原 OQ 质量值（Picard 默认为 true）。; `remove_duplicate_information` (可选): 清除重复标记。; `remove_alignment_information` (可选): 移除比对信息。; `sanitize` (可选): 移除无法一致还原的读段。; `sample_alias` (可选): 读组使用的样本名。; `library_name` (可选): 读组使用的文库名。; `sort_order` (可选): 输出排序方式。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为还原后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_revert_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_revert_sam(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 30) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_convert_sam_format()`

- **作用:** 把 SAM 转换为 BAM，或把 BAM 转换为 SAM（或 CRAM）。
- **处理方法:** 运行 Picard `gatk SamFormatConverter`。
- **API:** `gatk_convert_sam_format(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SamFormatConverter。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `output_path` (必需): 新建文件；格式由扩展名决定（.sam、.bam、.cram）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为转换后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_convert_sam_format

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_convert_sam_format(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 31) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_sam_to_fastq_with_tags()`

- **作用:** 把 SAM/BAM/CRAM 转换为 FASTQ，并同时写出由读段标签（例如 UMI）构成的 FASTQ。
- **处理方法:** 运行 Picard `gatk SamToFastqWithTags`。
- **API:** `gatk_sam_to_fastq_with_tags(input_path: 'str', fastq_path: 'str', sequence_tag_groups: 'Sequence[str]', *, quality_tag_groups: 'Sequence[str] | None' = None, second_end_fastq: 'str | None' = None, unpaired_fastq: 'str | None' = None, interleave: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SamToFastqWithTags。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM。; `fastq_path` (必需): 第一读段（交错时为全部读段）的新建 FASTQ。; `sequence_tag_groups` (必需): 其值写成额外 FASTQ 文件的标签；多个标签用逗号连接。; `quality_tag_groups` (可选): 与序列标签组配套的质量标签。; `second_end_fastq` (可选): 第二读段的新建 FASTQ。; `unpaired_fastq` (可选): 未配对读段的新建 FASTQ。; `interleave` (可选): 把一对读段写入同一个 FASTQ。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs` 列出所请求的 FASTQ 文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_sam_to_fastq_with_tags

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_sam_to_fastq_with_tags(
            input_path=str(input_path),
            fastq_path=str(folder / "fastq_path.out"),
            sequence_tag_groups=["RX"],
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 32) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_split_n_cigar_reads()`

- **作用:** 在 CIGAR 的 N 操作处把 RNA-seq 读段拆分为外显子片段并剪切悬挂部分。
- **处理方法:** 运行 `gatk SplitNCigarReads`。
- **API:** `gatk_split_n_cigar_reads(input_path: 'str', reference_path: 'str', output_path: 'str', *, max_bases_in_overhang: 'int | None' = None, max_mismatches_in_overhang: 'int | None' = None, skip_mapping_quality_transform: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SplitNCigarReads。
- **输入:** `input_path` (必需): CIGAR 含 N 操作的 RNA-seq BAM/CRAM/SAM。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的文件。; `max_bases_in_overhang` (可选): 拆分片段在被硬剪切前允许的最大悬挂碱基数。; `max_mismatches_in_overhang` (可选): 悬挂区允许的最大错配数。; `skip_mapping_quality_transform` (可选): 不把 MAPQ 255 改为 60。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为拆分后的文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_split_n_cigar_reads

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_split_n_cigar_reads(
            input_path=str(input_path),
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 33) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_split_sam()`

- **作用:** 把 SAM/BAM/CRAM 文件拆分为多个文件，按文库名或按读段数拆分。
- **处理方法:** 运行 Picard `SplitSamByLibrary`（`by=library`，无文库的读段写入 `unknown`）或 `SplitSamByNumberOfReads`（`by=number_of_reads`）；写出的文件列在 `files` 中。 `by`: `library` = SplitSamByLibrary，`number_of_reads` = SplitSamByNumberOfReads。
- **API:** `gatk_split_sam(input_path: 'str', output_dir: 'str', *, by: 'str' = 'library', split_to_n_files: 'int | None' = None, split_to_n_reads: 'int | None' = None, total_reads_in_input: 'int | None' = None, out_prefix: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SplitSamByLibrary, SplitSamByNumberOfReads。
- **输入:** `input_path` (必需): 输入 SAM/BAM/CRAM（`number_of_reads` 需按查询名分组）。; `output_dir` (必需): 存放拆分文件的新目录。; `split_to_n_files` (可选): 输出文件数。; `split_to_n_reads` (可选): 每个输出文件的读段数。; `total_reads_in_input` (可选): 读段总数，可省去一次计数遍历。; `out_prefix` (可选): 文件名前缀（默认 `shard`）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段、`mode` 和 `files`（`output_dir` 内的文件名）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_split_sam

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_split_sam(
            input_path=str(input_path),
            output_dir=str(folder / "output_dir"),
            by='library',
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 34) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_transfer_read_tags()`

- **作用:** 把未比对 SAM/BAM 中的读段标签（例如 UMI）复制到比对文件中对应的读段上。
- **处理方法:** 运行 `gatk TransferReadTags`（实验性）。
- **API:** `gatk_transfer_read_tags(input_path: 'str', unmapped_sam: 'str', output_path: 'str', read_tags: 'Sequence[str]', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：TransferReadTags。
- **输入:** `input_path` (必需): 接收标签的比对 BAM/SAM（按查询名排序）。; `unmapped_sam` (必需): 含标签的未比对 SAM/BAM（按查询名排序）。; `output_path` (必需): 新建的输出文件。; `read_tags` (必需): 要转移的标签名。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为新文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_transfer_read_tags

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    unmapped_sam = folder / "unmapped_sam.dat"; unmapped_sam.write_bytes(b"x")
    try:
        gatk_transfer_read_tags(
            input_path=str(input_path),
            unmapped_sam=str(unmapped_sam),
            output_path=str(folder / "output_path.out"),
            read_tags=["RX"],
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```

## 35) `dnakit.alignment.gatk_p01_alignment_cli_functions.gatk_unmark_duplicates()`

- **作用:** 清除所有读段的重复标记（0x400）。
- **处理方法:** 运行 `gatk UnmarkDuplicates`。
- **API:** `gatk_unmark_duplicates(input_path: 'str', output_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：UnmarkDuplicates。
- **输入:** `input_path` (必需): 输入 BAM/CRAM/SAM。; `output_path` (必需): 新建的文件。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为新文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.alignment.gatk_p01_alignment_cli_functions import gatk_unmark_duplicates

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_unmark_duplicates(
            input_path=str(input_path),
            output_path=str(folder / "output_path.out"),
            gatk_path="/nonexistent/gatk",
            verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出:**

```text
BackendUnavailableError GATK_UNAVAILABLE
```
