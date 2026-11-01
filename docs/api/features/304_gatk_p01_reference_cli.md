# GATK 参考序列、FASTA 与区间工具（GATK，P01 包）

这组 15 个函数位于 `dnakit.sequence.gatk_p01_reference_cli_functions`，调用 Broad Institute 的 [GATK](https://github.com/broadinstitute/gatk) 4.7.0.0 命令行（含其内置的 Picard 工具）。每个函数运行一个固定工具，或由模式参数选择的一族工具；不经过 shell，带超时和输出大小限制。`options` 只接受 `gatk <Tool> --help` 列出的长选项，已由参数设置的选项不能重复，GATK 会执行脚本的 `JAVASCRIPT_FILE` 不可用。输出写到显式的新路径，不会覆盖已有文件。

GATK 是外部程序：请安装（例如 `conda install -c bioconda gatk4` 或发布包）并准备 Java 17，然后把 `gatk` 启动器放进 PATH、设置 `DNAKIT_GATK` 或传入 `gatk_path`。找不到时抛出带安装提示的 `BackendUnavailableError`（代码 `GATK_UNAVAILABLE`）；工具运行失败抛出 `BackendExecutionError`（代码 `GATK_EXECUTION_FAILED`，附日志尾部）；工具返回 0 却未写出所请求的输出时抛出 `GATK_OUTPUT_MISSING`。DNAKit 不捆绑 GATK。Spark 工具以本地模式运行（可通过 `options` 设置 `spark-master`）。参考 FASTA 旁需要 `.fai` 和 `.dict`，BAM/VCF 需要各自索引，与 GATK 本身要求一致。结果为字典，每个函数的小节列出其额外字段。

## 1) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_index_feature_file()`

- **作用:** 为特征文件创建 tribble 索引（.idx；bgzip 压缩的 VCF 为 .tbi）。
- **处理方法:** 运行 `gatk IndexFeatureFile`；未给出 `index_path` 时，运行后检测输入旁新生成的索引并报告为 `outputs.index`。
- **API:** `gatk_index_feature_file(feature_path: 'str', *, index_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：IndexFeatureFile。
- **输入:** `feature_path` (必需): 特征文件（VCF、BED、GFF 等，tribble 支持的格式）。; `index_path` (可选): 显式指定索引文件；默认写在特征文件旁边。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.index` 为索引。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_index_feature_file

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    feature_path = folder / "feature_path.dat"; feature_path.write_bytes(b"x")
    try:
        gatk_index_feature_file(
            feature_path=str(feature_path),
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

## 2) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_sort_gff()`

- **作用:** 对 GFF3 文件排序并添加 flush（`###`）指令。
- **处理方法:** 运行 Picard `gatk SortGff`；给出 `sequence_dictionary` 时按其 contig 顺序排序。
- **API:** `gatk_sort_gff(input_path: 'str', output_path: 'str', *, sequence_dictionary: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：SortGff。
- **输入:** `input_path` (必需): GFF3 文件。; `output_path` (必需): 新建的已排序 GFF3 文件（含 flush 指令）。; `sequence_dictionary` (可选): 定义 contig 顺序的序列字典。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为排序后的 GFF3。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_sort_gff

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_sort_gff(
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

## 3) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_gtf_to_bed()`

- **作用:** 把 Gencode GTF 文件转换为 BED 文件。
- **处理方法:** 运行 `gatk GtfToBed`。
- **API:** `gatk_gtf_to_bed(gtf_path: 'str', output_path: 'str', sequence_dictionary: 'str', *, sort_by_transcript: 'bool | None' = None, use_basic_transcript: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：GtfToBed。
- **输入:** `gtf_path` (必需): Gencode GTF 文件。; `output_path` (必需): 新建的 BED 文件。; `sequence_dictionary` (必需): 参考序列的序列字典（.dict），用于确定 contig 顺序。; `sort_by_transcript` (可选): 按转录本排序输出。; `use_basic_transcript` (可选): 只保留标记为 `basic` 的转录本。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 BED 文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_gtf_to_bed

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    gtf_path = folder / "gtf_path.dat"; gtf_path.write_bytes(b"x")
    sequence_dictionary = folder / "sequence_dictionary.dat"; sequence_dictionary.write_bytes(b"x")
    try:
        gatk_gtf_to_bed(
            gtf_path=str(gtf_path),
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

## 4) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_bait_designer()`

- **作用:** 根据目标区间设计杂交捕获（外显子捕获）用的寡核苷酸探针。
- **处理方法:** 运行 Picard `gatk BaitDesigner`；写入 `output_dir` 的文件列在 `files` 中。
- **API:** `gatk_bait_designer(targets: 'str', reference_path: 'str', design_name: 'str', output_dir: 'str', *, bait_size: 'int | None' = None, bait_offset: 'int | None' = None, design_strategy: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：BaitDesigner。
- **输入:** `targets` (必需): 目标区域的区间列表。; `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `design_name` (必需): 探针设计名称。; `output_dir` (必需): 存放设计文件的新目录。; `bait_size` (可选): 每条探针长度。; `bait_offset` (可选): 相邻探针的间隔。; `design_strategy` (可选): 探针在目标上的布局策略。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `files`（`output_dir` 内的文件名）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_bait_designer

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    targets = folder / "targets.dat"; targets.write_bytes(b"x")
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_bait_designer(
            targets=str(targets),
            reference_path=str(reference_path),
            design_name="design1",
            output_dir=str(folder / "output_dir"),
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

## 5) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_bwa_mem_index_image()`

- **作用:** 为参考序列创建供 GATK BWA 工具使用的 BWA-MEM 索引镜像（`.img`）。
- **处理方法:** 运行 `gatk BwaMemIndexImageCreator`；镜像 `<fasta>.img` 写在 FASTA 旁边，且事先不能存在，报告为 `outputs.image`。
- **API:** `gatk_bwa_mem_index_image(reference_path: 'str', *, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：BwaMemIndexImageCreator。
- **输入:** `reference_path` (必需): 参考 FASTA。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.image` 为索引镜像。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_bwa_mem_index_image

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_bwa_mem_index_image(
            reference_path=str(reference_path),
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

## 6) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_check_reference_compatibility()`

- **作用:** 检查 BAM/CRAM 或 VCF 与一个或多个参考 FASTA 是否兼容。
- **处理方法:** 运行 `gatk CheckReferenceCompatibility`（实验性）；`alignment_path` 与 `variant_path` 至少给一个。
- **API:** `gatk_check_reference_compatibility(references_to_compare: 'Sequence[str]', *, alignment_path: 'str | None' = None, variant_path: 'str | None' = None, output_path: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CheckReferenceCompatibility。
- **输入:** `references_to_compare` (必需): 要测试的参考 FASTA 文件。; `alignment_path` (可选): 检查其序列字典的 BAM/CRAM。; `variant_path` (可选): 检查其 contig 头的 VCF。; `output_path` (可选): 新建的报告文件（默认写到标准输出）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段（未给 `output_path` 时表格在 `log_tail` 中）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_check_reference_compatibility

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    references_to_compare = folder / "references_to_compare.dat"; references_to_compare.write_bytes(b"x")
    alignment_path = folder / "alignment_path.dat"; alignment_path.write_bytes(b"x")
    try:
        gatk_check_reference_compatibility(
            references_to_compare=[str(references_to_compare)],
            alignment_path=str(alignment_path),
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

## 7) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_compare_references()`

- **作用:** 逐 contig 比较参考 FASTA 文件（MD5、长度、碱基差异）并汇总差异。
- **处理方法:** 运行 `gatk CompareReferences`（实验性）并解析其输出的表格。
- **API:** `gatk_compare_references(reference_path: 'str', references_to_compare: 'Sequence[str]', *, display_only_differing_sequences: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CompareReferences。
- **输入:** `reference_path` (必需): 主参考 FASTA。; `references_to_compare` (必需): 要比较的参考 FASTA 文件。; `display_only_differing_sequences` (可选): 只列出有差异的序列。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段，加上 `comparison_table`（打印表格的各行）。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_compare_references

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    references_to_compare = folder / "references_to_compare.dat"; references_to_compare.write_bytes(b"x")
    try:
        gatk_compare_references(
            reference_path=str(reference_path),
            references_to_compare=[str(references_to_compare)],
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

## 8) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_compose_str_table()`

- **作用:** 生成全基因组短串联重复（STR）位置表，用于 DRAGstr 模型校准。
- **处理方法:** 运行 `gatk ComposeSTRTableFile`。
- **API:** `gatk_compose_str_table(reference_path: 'str', output_path: 'str', *, max_period: 'int | None' = None, max_repeats: 'int | None' = None, decimation: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ComposeSTRTableFile。
- **输入:** `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的 STR 表（二进制）。; `max_period` (可选): 重复单元最大长度。; `max_repeats` (可选): 最大重复次数。; `decimation` (可选): 抽稀表（`period:factor,...`）。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 STR 表。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_compose_str_table

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_compose_str_table(
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

## 9) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_count_reference_bases()`

- **作用:** 统计 FASTA 文件的碱基：A、C、G、T、N 各自的数量，或非 N 碱基总数。
- **处理方法:** 运行 `gatk CountBasesInReference`（`mode=count_bases`，解析为 `base_counts`）或 Picard `NonNFastaSize`（`mode=non_n_size`，结果在输出文件中，可限定区间列表）。 `mode`: `count_bases` = CountBasesInReference，`non_n_size` = NonNFastaSize。
- **API:** `gatk_count_reference_bases(fasta_path: 'str', *, mode: 'str' = 'count_bases', output_path: 'str | None' = None, intervals: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：CountBasesInReference, NonNFastaSize。
- **输入:** `fasta_path` (必需): FASTA 文件。; `output_path` (部分模式必需): 非 N 碱基数的新建文件。; `intervals` (可选): 限制计数范围的区间列表。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`count_bases` 返回 `base_counts` 与 `total_bases`，`non_n_size` 返回 `outputs.output_path`。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_count_reference_bases

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    fasta_path = folder / "fasta_path.dat"; fasta_path.write_bytes(b"x")
    try:
        gatk_count_reference_bases(
            fasta_path=str(fasta_path),
            mode='count_bases',
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

## 10) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_extract_reference_sequences()`

- **作用:** 把参考序列中若干区间的序列提取到新的 FASTA 文件。
- **处理方法:** 运行 Picard `ExtractSequences`（`source=interval_list`，区间来自区间列表文件）或 `gatk FastaReferenceMaker`（`source=intervals`，区间写成 `chr:start-end`）。 `source`: `interval_list` = ExtractSequences，`intervals` = FastaReferenceMaker。
- **API:** `gatk_extract_reference_sequences(reference_path: 'str', output_path: 'str', *, source: 'str' = 'interval_list', interval_list: 'str | None' = None, intervals: 'Sequence[str] | None' = None, line_length: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ExtractSequences, FastaReferenceMaker。
- **输入:** `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的 FASTA 文件。; `interval_list` (部分模式必需): 要提取的区间列表（Picard 格式）。; `intervals` (部分模式必需): 要提取的区间，例如 `chr1:100-200`。; `line_length` (可选): 输出 FASTA 的行宽。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段和 `mode`；`outputs.output_path` 为 FASTA。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_extract_reference_sequences

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    interval_list = folder / "interval_list.dat"; interval_list.write_bytes(b"x")
    try:
        gatk_extract_reference_sequences(
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            interval_list=str(interval_list),
            source='interval_list',
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

## 11) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_fasta_alternate_reference_maker()`

- **作用:** 把 VCF 中的变异应用到参考 FASTA 上，生成替代参考序列。
- **处理方法:** 运行 `gatk FastaAlternateReferenceMaker`；可用 `intervals` 选项限制区域。
- **API:** `gatk_fasta_alternate_reference_maker(reference_path: 'str', variant_path: 'str', output_path: 'str', *, snp_mask: 'str | None' = None, use_iupac_sample: 'str | None' = None, line_width: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：FastaAlternateReferenceMaker。
- **输入:** `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `variant_path` (必需): 其等位基因将被应用的已索引 VCF。; `output_path` (必需): 新建的替代参考 FASTA。; `snp_mask` (可选): 要用 N 屏蔽的位置 VCF。; `use_iupac_sample` (可选): 杂合基因型用 IUPAC 码表示的样本。; `line_width` (可选): 输出 FASTA 的行宽。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 FASTA。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_fasta_alternate_reference_maker

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    variant_path = folder / "variant_path.dat"; variant_path.write_bytes(b"x")
    try:
        gatk_fasta_alternate_reference_maker(
            reference_path=str(reference_path),
            variant_path=str(variant_path),
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

## 12) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_find_bad_genomic_kmers_spark()`

- **作用:** 找出参考序列中高频出现的 k-mer（例如用于 SV 发现时屏蔽）。
- **处理方法:** 运行 `gatk FindBadGenomicKmersSpark`（beta，Spark 本地模式）。
- **API:** `gatk_find_bad_genomic_kmers_spark(reference_path: 'str', output_path: 'str', *, k_size: 'int | None' = None, kmer_max_dust_score: 'int | None' = None, high_copy_fasta: 'str | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：FindBadGenomicKmersSpark。
- **输入:** `reference_path` (必需): 参考 FASTA（也支持 .2bit）。; `output_path` (必需): 新建的高频 k-mer 文本文件。; `k_size` (可选): k-mer 长度。; `kmer_max_dust_score` (可选): 所报告 k-mer 的最大 DUST 分值。; `high_copy_fasta` (可选): 已知高拷贝序列的 FASTA。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为 k-mer 列表。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 此处未完整运行：参数和选项白名单已对照 4.7.0.0 帮助核对，但该工具需要的输入（GATK-SV 证据、Funcotator 数据源、三联体等）没有可用的玩具数据。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_find_bad_genomic_kmers_spark

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_find_bad_genomic_kmers_spark(
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

## 13) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_normalize_fasta()`

- **作用:** 重写 FASTA 文件，使所有序列行等长。
- **处理方法:** 运行 Picard `gatk NormalizeFasta`。
- **API:** `gatk_normalize_fasta(input_path: 'str', output_path: 'str', *, line_length: 'int | None' = None, truncate_sequence_names_at_whitespace: 'bool | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：NormalizeFasta。
- **输入:** `input_path` (必需): FASTA 文件。; `output_path` (必需): 新建的规整化 FASTA。; `line_length` (可选): 序列行长（默认 100）。; `truncate_sequence_names_at_whitespace` (可选): 序列名在第一个空白处截断。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为新 FASTA。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_normalize_fasta

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    input_path = folder / "input_path.dat"; input_path.write_bytes(b"x")
    try:
        gatk_normalize_fasta(
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

## 14) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_scatter_intervals_by_ns()`

- **作用:** 在参考序列的 N 连续段处切分，写出区间列表。
- **处理方法:** 运行 Picard `gatk ScatterIntervalsByNs`。
- **API:** `gatk_scatter_intervals_by_ns(reference_path: 'str', output_path: 'str', *, output_type: 'str | None' = None, max_to_merge: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ScatterIntervalsByNs。
- **输入:** `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的区间列表。; `output_type` (可选): 写出哪类区间。; `max_to_merge` (可选): 把不超过该长度的 N 段并入相邻 ACGT 区间。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs.output_path` 为区间列表。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_scatter_intervals_by_ns

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_scatter_intervals_by_ns(
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

## 15) `dnakit.sequence.gatk_p01_reference_cli_functions.gatk_shift_fasta()`

- **作用:** 生成移位后的 FASTA（环状 contig 如 chrM 旋转一个偏移量），并生成把坐标移回的文件。
- **处理方法:** 运行 `gatk ShiftFasta`（beta）。
- **API:** `gatk_shift_fasta(reference_path: 'str', output_path: 'str', shift_back_output: 'str', *, shift_offset_list: 'Sequence[str] | None' = None, interval_file_name: 'str | None' = None, line_width: 'int | None' = None, options: 'Mapping[str, Any] | None' = None, gatk_path: 'str | None' = None, java_max_memory: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **依赖:** 需要外部启动器 `gatk`（测试版本 GATK 4.7.0.0，Java 17），放在 PATH、`DNAKIT_GATK` 或通过 `gatk_path` 传入；生成本页的环境未安装 GATK，所以可运行示例展示已验证的“找不到启动器”行为（`BackendUnavailableError`，代码 `GATK_UNAVAILABLE`）。对应 GATK 工具：ShiftFasta。
- **输入:** `reference_path` (必需): 参考 FASTA（需 .fai 与 .dict）。; `output_path` (必需): 新建的移位 FASTA。; `shift_back_output` (必需): 把移位坐标映射回原坐标的新建 chain 文件。; `shift_offset_list` (可选): 每个 contig 的偏移量，顺序与字典一致。; `interval_file_name` (可选): 移位 contig 的新建区间文件。; `line_width` (可选): 输出 FASTA 的行宽。; `options`：该工具的 GATK 长选项 `{名称: 值}`（只接受 `gatk <Tool> --help` 列出的选项；`True`/`False` 转为 `true`/`false`，列表会重复该选项）；`gatk_path`、`java_max_memory`（如 `'4g'`）、`verbose`。.
- **输出:** 通用结果字段；`outputs` 列出 FASTA 和 shift-back 文件。 通用字段：`status`、`backend`、`operation`、`tool`、`version`、`command`、`outputs`、`log_tail`、`provenance`。
- **验证:** 已验证：在玩具数据上用 GATK 4.7.0.0 完整运行（见本包 verification 目录）。
- **示例代码:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.sequence.gatk_p01_reference_cli_functions import gatk_shift_fasta

with tempfile.TemporaryDirectory() as tmp:
    folder = Path(tmp)
    reference_path = folder / "reference_path.dat"; reference_path.write_bytes(b"x")
    try:
        gatk_shift_fasta(
            reference_path=str(reference_path),
            output_path=str(folder / "output_path.out"),
            shift_back_output=str(folder / "shift_back_output.out"),
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
