# SeqKit 与 screed：FASTA 和 FASTQ 命令行操作

这一组 `dnakit.fasta_fastq` 中的函数封装了两个处理序列文件的上游工具：外部命令行程序 [SeqKit](https://github.com/shenwei356/seqkit)（21 个函数，位于 `seqkit_cli_functions`）和 Python 包 [screed](https://github.com/dib-lab/screed)（3 个函数，位于 `screed_gap_functions`）。SeqKit 函数不经过 shell，只运行白名单中的一个子命令（`subseq`、`stats`、`sum`、`watch`、`fa2fq`、`fx2tab`、`fish`、`grep`、`duplicate`、`range`、`head-genome`、`split`、`concat`、`rename`、`replace`、`bam`、`merge-slides`），带超时和输出大小限制，并把结果写到明确指定的新路径；测试所用的 SeqKit 版本是 2.14.0。screed 函数把 FASTA 或 FASTQ 文件索引成 screed 数据库，并能把数据库导回文本。导入方式为 `from dnakit.fasta_fastq import ...`。

依赖都是可选的，只有调用函数时才会导入或查找。SeqKit 是外部程序：先安装（例如 `conda install -c bioconda seqkit` 或下载发布版二进制），再放入 PATH 或传入 `seqkit_path`。找不到时抛出 `BackendUnavailableError`（代码 `SEQKIT_UNAVAILABLE`）并给出安装提示；SeqKit 运行出错退出时抛出 `BackendExecutionError`（代码 `SEQKIT_EXECUTION_FAILED`）并附日志末尾。screed 需要额外依赖 `dnakit[bio-screed]`；其包元数据标注的许可证为 BSD-3-Clause。DNAKit 不内置 SeqKit。

每个 SeqKit 函数都会先校验参数，绝不覆盖已有输出，并返回字典。所有 SeqKit 结果共有这些字段：`status`（`completed`）、`backend`（`seqkit`）、`operation`、`version`、`command`、`outputs`、`log_tail` 和 `provenance`；各函数小节只列出额外字段。生成本页所用的环境没有安装 SeqKit，因此 21 个 SeqKit 示例展示的是在任何机器上都真实成立的参数校验或缺少可执行文件时的行为，每个示例后面都附有安装 SeqKit 后运行的调用。三个 screed 示例是真实运行的。

## 示例数据

SeqKit 示例在临时目录中自己创建很小的输入。screed 示例读取下面两个文件。

**示例文件：** `reads.fa`

```text
>r1 first read
ACGTACGT
>r2 second read
GGCCTTAA
```

**示例文件：** `reads.fq`

```text
@r1 first read
ACGTACGT
+
IIIIHHHH
@r2 second read
GGCCTTAA
+
########
```

## 1) `dnakit.fasta_fastq.seqkit_bam_alignment_stats()`

- **作用：** 统计 BAM 文件中比对上、未比对、主比对、次要比对和补充比对的数量。
- **处理方法：** 运行 `seqkit bam --stat`（全文件扫描），或在 `index_based` 为真时运行 `seqkit bam --idx-stat`，并加 `--quiet-mode`，然后解析 SeqKit 输出到 stderr 的制表符分隔表格；扫描模式的第一列是 `PrimAlnPerc`，索引模式是 `AlnPerc`。
- **API：** `seqkit_bam_alignment_stats(bam_path: 'str', *, index_based: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：不存在的 `seqkit_path` 会抛出 `BackendUnavailableError`（代码 `SEQKIT_UNAVAILABLE`）。
- **输入：** `bam_path` 是已存在的 BAM 文件。`index_based` 必须是布尔值；为 True 时读取 BAM 索引，需要 BAM 旁有 `.bai` 文件。
- **输出：** 通用结果字典，另有 `records`（解析后的表格行，数值已转换，`File` 列保留为文本）和 `index_based`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.fasta_fastq import seqkit_bam_alignment_stats

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    try:
        seqkit_bam_alignment_stats(str(bam), index_based=True, seqkit_path="/nonexistent/seqkit", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
```

- **示例输出：**

```text
BackendUnavailableError SEQKIT_UNAVAILABLE
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_bam_alignment_stats

result = seqkit_bam_alignment_stats("reads.bam", verbose=False)
print(result["records"])
```

## 2) `dnakit.fasta_fastq.seqkit_bam_bundles()`

- **作用：** 把 BAM 文件按位点或按不少于指定读段数的捆绑（bundle）切分。
- **处理方法：** 先创建 `output_dir`，再以该目录为工作目录运行 `seqkit bam --bundle=<min_reads> --quiet-mode`，解析 SeqKit 日志中的 `Bundle` 表格以及 “Written N BAM records to L loci and B bundles” 汇总行。SeqKit 把各捆绑写入 `output_dir` 内的 `<bam 文件名>_bundles/`。
- **API：** `seqkit_bam_bundles(bam_path: 'str', output_dir: 'str', *, min_reads: 'int' = -1, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：不存在的 `seqkit_path` 会抛出 `BackendUnavailableError`（代码 `SEQKIT_UNAVAILABLE`），而此时输出目录已经创建。
- **输入：** `bam_path` 是已存在的 BAM 文件。`output_dir` 必须尚不存在，其父目录必须存在。`min_reads` 是不小于 -1 的整数；默认的 -1 表示每个位点一个捆绑。
- **输出：** 通用结果字典，另有 `records`（捆绑表格行，`Chrom` 保留为文本）、`min_reads`，以及在找到汇总行时的 `record_count`、`locus_count` 和 `bundle_count`。`outputs` 包含 `directory` 和 `bundles`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.fasta_fastq import seqkit_bam_bundles

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    out = Path(tmp) / "bundles"
    try:
        seqkit_bam_bundles(str(bam), str(out), min_reads=10, seqkit_path="/nonexistent/seqkit", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
    print("output directory created before the failure:", out.is_dir())
    try:
        seqkit_bam_bundles(str(bam), str(out), seqkit_path="/nonexistent/seqkit", verbose=False)
    except Exception as error:
        print(type(error).__name__, error.code)
```

- **示例输出：**

```text
BackendUnavailableError SEQKIT_UNAVAILABLE
output directory created before the failure: True
ConfigurationError CONFIGURATION_ERROR
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_bam_bundles

result = seqkit_bam_bundles("reads.bam", "bundles_out", min_reads=100, verbose=False)
print(result["record_count"], result["locus_count"], result["bundle_count"])
```

## 3) `dnakit.fasta_fastq.seqkit_bam_feature_histogram()`

- **作用：** 给 BAM 中每条记录的某个特征（如 `ReadLen`、`MapQual`、`MeanQual`）生成直方图。
- **处理方法：** 运行 `seqkit bam --field=<field> --bins=<bins> --dump --quiet-mode`（按需加 `--log` 和 `--img=<image_path>`），并解析 SeqKit 打印的 `Bin` 表格。函数只检查字段名的形状，是否为已知特征由 SeqKit 判断（`seqkit bam -H` 可列出）。按源码注释，SeqKit 2.14.0 在 BAM 没有 NM 标签时会在 `Acc` 上崩溃，所以不传 `Acc`。
- **API：** `seqkit_bam_feature_histogram(bam_path: 'str', field: 'str', *, bins: 'int' = -1, log_transform: 'bool' = False, image_path: 'str | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `bam_path` 是已存在的 BAM 文件。`field` 由字母和数字组成且以字母开头。`bins` 是不小于 -1 的整数（-1 由 SeqKit 自行决定）。`log_transform` 是布尔值。`image_path` 可选，必须是已存在文件夹中的新文件。
- **输出：** 通用结果字典，另有 `records`（直方图各区间）和 `field`。请求了图像时，`outputs` 含 `image`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_bam_feature_histogram

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    try:
        seqkit_bam_feature_histogram(str(bam), "Read Len", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] field must be a SeqKit BAM feature name such as ReadLen.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_bam_feature_histogram

result = seqkit_bam_feature_histogram("reads.bam", "ReadLen", bins=10, verbose=False)
print(result["field"], len(result["records"]))
```

## 4) `dnakit.fasta_fastq.seqkit_bam_top_records()`

- **作用：** 把某个特征值最大的若干条 BAM 记录保存到新的 BAM 文件。
- **处理方法：** 运行 `seqkit bam --field=<field> --top-bam=<output_bam> --top-size=<top_size> --quiet-mode`，并在结束后检查新 BAM 是否存在。
- **API：** `seqkit_bam_top_records(bam_path: 'str', output_bam: 'str', *, field: 'str', top_size: 'int' = 100, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `bam_path` 是已存在的 BAM 文件。`output_bam` 必须不存在，其所在文件夹必须存在。`field` 是特征名（字母和数字，以字母开头），必填。`top_size` 是不小于 1 的整数（默认 100）。
- **输出：** 通用结果字典，另有 `field` 和 `top_size`。`outputs` 含 `bam`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_bam_top_records

with tempfile.TemporaryDirectory() as tmp:
    bam = Path(tmp) / "reads.bam"
    bam.write_bytes(b"not a real BAM")
    try:
        seqkit_bam_top_records(str(bam), str(Path(tmp) / "top.bam"), field="ReadLen", top_size=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] top_size must be at least 1.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_bam_top_records

result = seqkit_bam_top_records("reads.bam", "longest.bam", field="ReadLen", top_size=10, verbose=False)
print(result["outputs"]["bam"].endswith("longest.bam"))
```

## 5) `dnakit.fasta_fastq.seqkit_concat_by_id()`

- **作用：** 把多个文件中 ID 相同的序列按文件顺序首尾连接。
- **处理方法：** 运行 `seqkit concat`，带 `--separator=<separator>`，`full_join` 为真时再加 `--full`（以及 `--fill=<fill>`），用 `-o` 写入 `output_path`，之后统计输出的记录数。
- **API：** `seqkit_concat_by_id(input_paths: 'Sequence[str]', output_path: 'str', *, full_join: 'bool' = False, fill: 'str | None' = None, separator: 'str' = '|', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_paths` 是非空的已存在文件列表。`output_path` 必须不存在。`full_join` 是布尔值。`fill` 可选，只能与 `full_join=True` 同用，且只能包含字母或 `-`。`separator` 是连接描述文字的字符串（默认 `|`）。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_concat_by_id

with tempfile.TemporaryDirectory() as tmp:
    a = Path(tmp) / "a.fa"
    b = Path(tmp) / "b.fa"
    a.write_text(">s1\nAAAA\n")
    b.write_text(">s1\nCCCC\n")
    try:
        seqkit_concat_by_id([str(a), str(b)], str(Path(tmp) / "joined.fa"), fill="N", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] fill only applies with full_join=True.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_concat_by_id

result = seqkit_concat_by_id(["part1.fa", "part2.fa"], "joined.fa", full_join=True, fill="N", verbose=False)
print(result["record_count"])
```

## 6) `dnakit.fasta_fastq.seqkit_count_matching_records()`

- **作用：** 统计 ID、完整名称或序列与任一模式匹配的记录数。
- **处理方法：** 运行 `seqkit grep --count`，每个模式对应一个 `--pattern=`，并按需加 `--by-name`、`--by-seq`、`--use-regexp`、`--invert-match` 和 `--ignore-case`，把输出中最后一个纯数字行读作计数。不写任何文件。
- **API：** `seqkit_count_matching_records(input_path: 'str', *, patterns: 'Sequence[str]', by_name: 'bool' = False, by_sequence: 'bool' = False, use_regexp: 'bool' = False, invert_match: 'bool' = False, ignore_case: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`patterns` 是非空的、无逗号的非空字符串列表（SeqKit 以逗号拆分模式列表）。其余选项是布尔值。
- **输出：** 通用结果字典，另有 `count`（整数）。`outputs` 为空。若 SeqKit 没有输出计数，则抛出代码为 `SEQKIT_OUTPUT_UNREADABLE` 的 `BackendExecutionError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_count_matching_records

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n>s2\nCCCC\n")
    try:
        seqkit_count_matching_records(str(fasta), patterns=["s1,s2"], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] each pattern must be a non-empty string without commas (SeqKit splits pattern lists on commas).
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_count_matching_records

result = seqkit_count_matching_records("reads.fa", patterns=["r1", "r2"], verbose=False)
print(result["count"])
```

## 7) `dnakit.fasta_fastq.seqkit_duplicate_records()`

- **作用：** 把每条记录原地重复写出若干次，ID 保持不变。
- **处理方法：** 运行 `seqkit duplicate --times=<times>`，用 `-o` 写入 `output_path`，之后统计输出的记录数。
- **API：** `seqkit_duplicate_records(input_path: 'str', output_path: 'str', *, times: 'int', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`output_path` 必须不存在。`times` 是必填的、不小于 1 的整数。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_duplicate_records

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n")
    try:
        seqkit_duplicate_records(str(fasta), str(Path(tmp) / "dup.fa"), times=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] times must be at least 1.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_duplicate_records

result = seqkit_duplicate_records("reads.fa", "tripled.fa", times=3, verbose=False)
print(result["record_count"])
```

## 8) `dnakit.fasta_fastq.seqkit_extended_stats()`

- **作用：** 为每个文件计算扩展统计：长度四分位数、N50 与 Nx、Q20/Q30、平均质量和 GC。
- **处理方法：** 运行 `seqkit stats --all --tabular --fq-encoding=<fastq_encoding>`（`n_values` 的每个值再加一个 `--N=<value>`），把表格写入 `output_path`，再读回该文件，以第一行为表头并转换数值单元格。
- **API：** `seqkit_extended_stats(input_paths: 'Sequence[str]', output_path: 'str', *, fastq_encoding: 'str' = 'sanger', n_values: 'Sequence[int]' = (), seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_paths` 是非空的已存在文件列表。`output_path` 必须不存在。`fastq_encoding` 取 `sanger`、`solexa`、`illumina-1.3+`、`illumina-1.5+`、`illumina-1.8+` 之一。`n_values` 是 0 到 100 的整数序列。
- **输出：** 通用结果字典，另有 `records`，每个输入文件一个字典，含各表格列（`file`、`format`、`type` 保留为文本）。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_extended_stats

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    try:
        seqkit_extended_stats([str(fastq)], str(Path(tmp) / "stats.tsv"), fastq_encoding="phred64", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
    try:
        seqkit_extended_stats([str(fastq)], str(Path(tmp) / "stats.tsv"), n_values=[101], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] fastq_encoding must be one of ['sanger', 'solexa', 'illumina-1.3+', 'illumina-1.5+', 'illumina-1.8+']; got 'phred64'.
ConfigurationError [CONFIGURATION_ERROR] n_values must be at most 100.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_extended_stats

result = seqkit_extended_stats(["reads.fq"], "stats.tsv", n_values=[90], verbose=False)
print(result["records"][0]["num_seqs"])
```

## 9) `dnakit.fasta_fastq.seqkit_fish_search()`

- **作用：** 用局部比对在较长的目标序列中定位较短的查询序列。
- **处理方法：** 对目标文件运行 `seqkit fish --query-fastx=<query_path> --min-qual=<min_mapping_quality>`（按需加 `--aln-params`、`--all`、`--stranded`）。由于 SeqKit 2.14.0 的 `fish` 忽略 `-o`，函数从 SeqKit 输出中取出以 `Ref` 为表头的表格，自行写入 `output_path`，再解析各行（`Ref`、`Query`、`Strand` 保留为文本）。
- **API：** `seqkit_fish_search(target_path: 'str', query_path: 'str', output_path: 'str', *, min_mapping_quality: 'float' = 5.0, search_all: 'bool' = False, stranded: 'bool' = False, alignment_params: 'str | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `target_path` 和 `query_path` 是已存在的 FASTA/FASTQ 文件。`output_path` 必须不存在。`min_mapping_quality` 是数值（默认 5.0）。`alignment_params` 为 None，或形如 `match,mismatch,gap_open,gap_extend` 的四个整数（SeqKit 默认 4,-4,-2,-1）。`search_all` 和 `stranded` 是布尔值。按源码文档字符串，目标含 N 时 SeqKit 2.14.0 以退出码 2 结束，会以 `BackendExecutionError` 报出。
- **输出：** 通用结果字典，另有 `records`（每个命中一行）。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_fish_search

with tempfile.TemporaryDirectory() as tmp:
    target = Path(tmp) / "target.fa"
    query = Path(tmp) / "query.fa"
    target.write_text(">t1\nACGTACGTACGTAAGGCCTT\n")
    query.write_text(">q1\nAAGGCC\n")
    try:
        seqkit_fish_search(str(target), str(query), str(Path(tmp) / "hits.tsv"), alignment_params="4,-4,-2", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] alignment_params must look like '4,-4,-2,-1'.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_fish_search

result = seqkit_fish_search("targets.fa", "queries.fa", "hits.tsv", search_all=True, verbose=False)
print(len(result["records"]))
```

## 10) `dnakit.fasta_fastq.seqkit_fastq_records_for_fasta()`

- **作用：** 返回在 FASTA 中有相同 ID 且序列相同的 FASTQ 记录。
- **处理方法：** 对 FASTQ 文件运行 `seqkit fa2fq --fasta-file=<fasta_path>`（按需加 `--only-positive-strand`），用 `-o` 写入 `output_path`。按源码文档字符串，SeqKit 2.14.0 默认在两条链上比较完整序列，除非 `only_positive_strand` 为 True；FASTA 中不存在或序列不同的 FASTQ 记录不会返回。
- **API：** `seqkit_fastq_records_for_fasta(fastq_path: 'str', fasta_path: 'str', output_path: 'str', *, only_positive_strand: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：不存在的 `seqkit_path` 会抛出 `BackendUnavailableError`（代码 `SEQKIT_UNAVAILABLE`），且不会写出输出文件。
- **输入：** `fastq_path` 和 `fasta_path` 是已存在的文件。`output_path` 必须不存在。`only_positive_strand` 是布尔值。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.fasta_fastq import seqkit_fastq_records_for_fasta

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fasta = Path(tmp) / "reads.fa"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    fasta.write_text(">r1\nACGT\n")
    out = Path(tmp) / "matched.fq"
    try:
        seqkit_fastq_records_for_fasta(str(fastq), str(fasta), str(out), seqkit_path="/nonexistent/seqkit", verbose=False)
    except BackendUnavailableError as error:
        print(type(error).__name__, error.code)
    print("output file written:", out.exists())
```

- **示例输出：**

```text
BackendUnavailableError SEQKIT_UNAVAILABLE
output file written: False
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_fastq_records_for_fasta

result = seqkit_fastq_records_for_fasta("reads.fq", "reads.fa", "matched.fq", verbose=False)
print(result["record_count"])
```

## 11) `dnakit.fasta_fastq.seqkit_flanking_sequences()`

- **作用：** 提取 BED 或 GTF 特征及其上游、下游侧翼序列。
- **处理方法：** 运行 `seqkit subseq`，带 `--bed=<feature_path>` 或 `--gtf=<feature_path>`，按需加 `--feature=<name>`（每个 GTF 特征类型一个）、`--up-stream`、`--down-stream` 和 `--only-flank`，用 `-o` 写入 `output_path`。SeqKit 可能在 FASTA 旁留下 `<序列文件>.seqkit.fai` 索引；若存在，其路径作为 `index_file` 返回。
- **API：** `seqkit_flanking_sequences(sequence_path: 'str', feature_path: 'str', output_path: 'str', *, feature_format: 'str' = 'bed', up_stream: 'int' = 0, down_stream: 'int' = 0, only_flank: 'bool' = False, feature_types: 'Sequence[str] | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `sequence_path` 和 `feature_path` 是已存在的文件。`output_path` 必须不存在。`feature_format` 为 `bed` 或 `gtf`。`up_stream` 和 `down_stream` 是不小于 0 的整数。`only_flank=True` 时二者中必须恰有一个大于 0。`feature_types`（名称列表）只允许用于 `gtf`。
- **输出：** 通用结果字典，另有 `record_count` 和 `index_file`（路径字符串或 None）。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_flanking_sequences

with tempfile.TemporaryDirectory() as tmp:
    genome = Path(tmp) / "genome.fa"
    bed = Path(tmp) / "genes.bed"
    genome.write_text(">chr1\nACGTACGTACGTACGTACGT\n")
    bed.write_text("chr1\t5\t10\tgeneA\n")
    try:
        seqkit_flanking_sequences(str(genome), str(bed), str(Path(tmp) / "out.fa"), up_stream=2, down_stream=2, only_flank=True, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
    try:
        seqkit_flanking_sequences(str(genome), str(bed), str(Path(tmp) / "out.fa"), feature_types=["gene"], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] only_flank needs exactly one of up_stream or down_stream above 0.
ConfigurationError [CONFIGURATION_ERROR] feature_types only applies to feature_format='gtf'.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_flanking_sequences

result = seqkit_flanking_sequences("genome.fa", "genes.bed", "genes_flanked.fa", up_stream=100, down_stream=100, verbose=False)
print(result["record_count"])
```

## 12) `dnakit.fasta_fastq.seqkit_grep_motif()`

- **作用：** 保留序列中含有任一给定基序的记录。
- **处理方法：** 运行 `seqkit grep --by-seq`，每个基序对应一个 `--pattern=`，并按需加 `--degenerate`、`--max-mismatch=<n>`、`--only-positive-strand`、`--ignore-case` 和 `--invert-match`，用 `-o` 写入 `output_path`。默认按字面匹配，并同时检查两条链，除非 `only_positive_strand` 为 True。
- **API：** `seqkit_grep_motif(input_path: 'str', output_path: 'str', *, motifs: 'Sequence[str]', degenerate: 'bool' = False, max_mismatch: 'int' = 0, only_positive_strand: 'bool' = False, ignore_case: 'bool' = False, invert_match: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`output_path` 必须不存在。`motifs` 是非空的、无逗号的非空字符串列表。`max_mismatch` 是不小于 0 的整数。其余选项是布尔值。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_grep_motif

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAGGTCTCAA\n>s2\nCCCCCCCC\n")
    try:
        seqkit_grep_motif(str(fasta), str(Path(tmp) / "hits.fa"), motifs=["GGTCTC,GAGACC"], verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] each pattern must be a non-empty string without commas (SeqKit splits pattern lists on commas).
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_grep_motif

result = seqkit_grep_motif("reads.fa", "with_site.fa", motifs=["GGTCTC"], verbose=False)
print(result["record_count"])
```

## 13) `dnakit.fasta_fastq.seqkit_head_genome()`

- **作用：** 返回文件中第一个基因组的记录；基因组指名称前导词相同的一段连续记录。
- **处理方法：** 运行 `seqkit head-genome --mini-common-words=<minimal_common_words>`，用 `-o` 写入 `output_path`，之后统计输出的记录数。
- **API：** `seqkit_head_genome(input_path: 'str', output_path: 'str', *, minimal_common_words: 'int' = 1, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`output_path` 必须不存在。`minimal_common_words` 是不小于 1 的整数（默认 1）。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_head_genome

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "genomes.fa"
    fasta.write_text(">Escherichia coli contig1\nACGT\n>Escherichia coli contig2\nGGCC\n")
    try:
        seqkit_head_genome(str(fasta), str(Path(tmp) / "first.fa"), minimal_common_words=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] minimal_common_words must be at least 1.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_head_genome

result = seqkit_head_genome("genomes.fa", "first_genome.fa", minimal_common_words=2, verbose=False)
print(result["record_count"])
```

## 14) `dnakit.fasta_fastq.seqkit_mean_read_quality()`

- **作用：** 报告每条 FASTQ 读段的平均质量。
- **处理方法：** 运行 `seqkit fx2tab --name --avg-qual --qual-ascii-base=<quality_ascii_base>`，用 `-o` 把两列表格写入 `output_path`，再解析成行。数值是 SeqKit 自己计算的平均值。
- **API：** `seqkit_mean_read_quality(fastq_path: 'str', output_path: 'str', *, quality_ascii_base: 'int' = 33, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `fastq_path` 是已存在的 FASTQ 文件。`output_path` 必须不存在。`quality_ascii_base` 是 0 到 127 的整数（默认 33，即 Phred+33）。
- **输出：** 通用结果字典，另有 `records`，每条读段一行 `{"id", "avg_qual"}`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_mean_read_quality

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    try:
        seqkit_mean_read_quality(str(fastq), str(Path(tmp) / "quality.tsv"), quality_ascii_base=200, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] quality_ascii_base must be at most 127.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_mean_read_quality

result = seqkit_mean_read_quality("reads.fq", "quality.tsv", verbose=False)
print(result["records"])
```

## 15) `dnakit.fasta_fastq.seqkit_merge_slides_to_bed()`

- **作用：** 把相互重叠的滑动窗口 ID 合并成 BED 区间。
- **处理方法：** 对窗口 ID 文本文件（如每行一个 `chr1_sliding:1-1000`）运行 `seqkit merge-slides --min-overlap=<min_overlap> --max-gap=<max_gap>`，用 `-o` 把 BED 输出写入 `output_path`，再解析成行。
- **API：** `seqkit_merge_slides_to_bed(window_names_path: 'str', output_path: 'str', *, min_overlap: 'int' = 1, max_gap: 'int' = 0, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `window_names_path` 是已存在的文本文件。`output_path` 必须不存在。`min_overlap` 是不小于 1 的整数（默认 1）。`max_gap` 是不小于 0 的整数（默认 0，表示不限制窗口起点之间的距离）。
- **输出：** 通用结果字典，另有 `records`，每个合并后的区间一行 `{"chrom", "start", "end"}`（起点 0 起始）。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_merge_slides_to_bed

with tempfile.TemporaryDirectory() as tmp:
    names = Path(tmp) / "windows.txt"
    names.write_text("chr1_sliding:1-10\nchr1_sliding:6-15\n")
    try:
        seqkit_merge_slides_to_bed(str(names), str(Path(tmp) / "merged.bed"), min_overlap=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] min_overlap must be at least 1.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_merge_slides_to_bed

result = seqkit_merge_slides_to_bed("windows.txt", "merged.bed", verbose=False)
print(result["records"])
```

## 16) `dnakit.fasta_fastq.seqkit_range_records()`

- **作用：** 按位置选取记录，位置可从文件开头或结尾计数。
- **处理方法：** 运行 `seqkit range --range=<start>:<end>`，用 `-o` 写入 `output_path`。位置从 1 开始且含两端；负数从末尾计数，所以 (1, 12) 是前 12 条记录，(-3, -1) 是最后三条。
- **API：** `seqkit_range_records(input_path: 'str', output_path: 'str', *, start: 'int', end: 'int', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`output_path` 必须不存在。`start` 和 `end` 是必填整数（布尔值会被拒绝）。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_range_records

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n>s2\nCCCC\n>s3\nGGGG\n")
    try:
        seqkit_range_records(str(fasta), str(Path(tmp) / "range.fa"), start="1", end=2, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] start must be an integer.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_range_records

result = seqkit_range_records("reads.fa", "last_two.fa", start=-2, end=-1, verbose=False)
print(result["record_count"])
```

## 17) `dnakit.fasta_fastq.seqkit_read_feature_histogram()`

- **作用：** 给 FASTA/FASTQ 文件的某个读段特征（长度、平均质量、GC 或 GC 偏斜）生成直方图。
- **处理方法：** 运行 `seqkit watch --fields=<field> --bins=<bins> --qual-ascii-base=<quality_ascii_base> --dump --quiet-mode`（按需加 `--log`、`--validate-seq` 和 `--img=<image_path>`），并解析 SeqKit 输出中的 `Bin` 表格。报告只在文件读完时生成一次，不提供在线监控选项。
- **API：** `seqkit_read_feature_histogram(input_path: 'str', field: 'str', *, bins: 'int' = -1, log_transform: 'bool' = False, quality_ascii_base: 'int' = 33, validate_seq: 'bool' = False, image_path: 'str | None' = None, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`field` 取 `ReadLen`、`MeanQual`、`GC`、`GCSkew` 之一。`bins` 是不小于 -1 的整数。`quality_ascii_base` 是 0 到 127 的整数。`log_transform` 和 `validate_seq` 是布尔值。`image_path` 可选，必须是已存在文件夹中的新文件。
- **输出：** 通用结果字典，另有 `records`（直方图各区间）和 `field`。请求了图像时，`outputs` 含 `image`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_read_feature_histogram

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGT\n+\nIIII\n")
    try:
        seqkit_read_feature_histogram(str(fastq), "Length", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] field must be one of ['ReadLen', 'MeanQual', 'GC', 'GCSkew']; got 'Length'.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_read_feature_histogram

result = seqkit_read_feature_histogram("reads.fq", "MeanQual", bins=5, verbose=False)
print(result["field"], len(result["records"]))
```

## 18) `dnakit.fasta_fastq.seqkit_rename_duplicate_ids()`

- **作用：** 给重复的 ID 加计数后缀，使其唯一。
- **处理方法：** 运行 `seqkit rename --start-num=<start_number> --separator=<separator>`（按需加 `--by-name` 和 `--rename-1st-rec`），用 `-o` 写入 `output_path`。除非 `rename_first` 为真，首次出现的 ID 保持不变，所以默认情况下 `id, id, id` 变成 `id, id_2, id_3`。
- **API：** `seqkit_rename_duplicate_ids(input_path: 'str', output_path: 'str', *, by_name: 'bool' = False, start_number: 'int' = 2, separator: 'str' = '_', rename_first: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`output_path` 必须不存在。`by_name` 和 `rename_first` 是布尔值。`start_number` 是不小于 1 的整数（默认 2）。`separator` 是非空字符串（默认 `_`）。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_rename_duplicate_ids

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "dups.fa"
    fasta.write_text(">s1\nAAAA\n>s1\nCCCC\n")
    try:
        seqkit_rename_duplicate_ids(str(fasta), str(Path(tmp) / "unique.fa"), start_number=0, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] start_number must be at least 1.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_rename_duplicate_ids

result = seqkit_rename_duplicate_ids("dups.fa", "unique.fa", verbose=False)
print(result["record_count"])
```

## 19) `dnakit.fasta_fastq.seqkit_replace_sequence_regex()`

- **作用：** 对文件中的序列做正则表达式替换。
- **处理方法：** 运行 `seqkit replace --by-seq --pattern=<pattern> --replacement=<replacement>`（按需加 `--ignore-case`），用 `-o` 写入 `output_path`。替换内容可以使用 `$1` 之类的捕获组。对 FASTQ，序列长度不能变化，因为质量字符串长度保持不变。
- **API：** `seqkit_replace_sequence_regex(input_path: 'str', output_path: 'str', *, pattern: 'str', replacement: 'str', ignore_case: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`output_path` 必须不存在。`pattern` 是非空字符串。`replacement` 是字符串（可为空）。`ignore_case` 是布尔值。
- **输出：** 通用结果字典，另有 `record_count`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_replace_sequence_regex

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAANNAA\n")
    try:
        seqkit_replace_sequence_regex(str(fasta), str(Path(tmp) / "fixed.fa"), pattern="N", replacement=None, verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] replacement must be a string.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_replace_sequence_regex

result = seqkit_replace_sequence_regex("reads.fa", "no_n.fa", pattern="N", replacement="A", verbose=False)
print(result["record_count"])
```

## 20) `dnakit.fasta_fastq.seqkit_sequence_digest()`

- **作用：** 为每个文件的序列计算与顺序相关的消息摘要。
- **处理方法：** 运行 `seqkit sum`（按需加 `--include-id`、`--remove-gaps` 和 `--rna2dna`），用 `-o` 把输出写入 `output_path`，并把每行解析为摘要和文件名。摘要形如 `seqkit.v0.1_...`，不是文件的 SHA-256。
- **API：** `seqkit_sequence_digest(input_paths: 'Sequence[str]', output_path: 'str', *, include_ids: 'bool' = False, remove_gaps: 'bool' = False, rna_to_dna: 'bool' = False, seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：已存在的输出文件会被拒绝且保持不变。
- **输入：** `input_paths` 是非空的已存在文件列表。`output_path` 必须不存在（已有文件绝不会被覆盖）。`include_ids`、`remove_gaps` 和 `rna_to_dna` 是布尔值。
- **输出：** 通用结果字典，另有 `records`，每个输入文件一行 `{"digest", "file"}`。`outputs` 含 `output`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_sequence_digest

with tempfile.TemporaryDirectory() as tmp:
    fasta = Path(tmp) / "reads.fa"
    fasta.write_text(">s1\nAAAA\n")
    existing = Path(tmp) / "digest.tsv"
    existing.write_text("old\n")
    try:
        seqkit_sequence_digest([str(fasta)], str(existing), verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error.code)
    print("existing file kept:", existing.read_text().strip())
```

- **示例输出：**

```text
ConfigurationError CONFIGURATION_ERROR
existing file kept: old
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_sequence_digest

result = seqkit_sequence_digest(["reads.fa"], "digest.tsv", verbose=False)
print(result["records"][0]["digest"].startswith("seqkit.v0.1_"))
```

## 21) `dnakit.fasta_fastq.seqkit_split_by_region()`

- **作用：** 按某个区域内的子序列（例如条形码）把记录拆分到不同文件。
- **处理方法：** 运行 `seqkit split --by-region=<region> --out-dir=<output_dir>`；SeqKit 创建该目录，并为每种不同的子序列写一个文件。之后函数列出目录中的文件并统计每个文件的记录数。
- **API：** `seqkit_split_by_region(input_path: 'str', output_dir: 'str', *, region: 'str', seqkit_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `seqkit`（测试版本为 2.14.0）位于 PATH 中，或显式传入 `seqkit_path`；生成本页的环境没有安装它，也不需要任何 Python 包。因此下面的示例展示不依赖 SeqKit 的行为：非法参数会在查找 SeqKit 之前抛出 `ConfigurationError`。
- **输入：** `input_path` 是已存在的 FASTA/FASTQ 文件。`output_dir` 必须不存在，其父目录必须存在。`region` 是形如 `start:end` 的字符串，可带负号，例如 `1:5` 或 `-12:-1`。
- **输出：** 通用结果字典，另有 `records`，每个写出的文件一行 `{"file", "record_count"}`。`outputs` 含 `directory`。若 SeqKit 没有创建该目录，则抛出代码为 `SEQKIT_OUTPUT_MISSING` 的 `BackendExecutionError`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import seqkit_split_by_region

with tempfile.TemporaryDirectory() as tmp:
    fastq = Path(tmp) / "reads.fq"
    fastq.write_text("@r1\nACGTAC\n+\nIIIIII\n")
    try:
        seqkit_split_by_region(str(fastq), str(Path(tmp) / "split"), region="first5", verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, error)
```

- **示例输出：**

```text
ConfigurationError [CONFIGURATION_ERROR] region must look like '1:5' or '-12:-1'.
```

**安装依赖后的用法：**

```python
from dnakit.fasta_fastq import seqkit_split_by_region

result = seqkit_split_by_region("reads.fq", "by_barcode", region="1:4", verbose=False)
print(result["records"])
```

## 22) `dnakit.fasta_fastq.create_screed_database()`

- **作用：** 把 FASTA 或 FASTQ 文件索引成 screed 数据库，路径由调用者指定。
- **处理方法：** 用 `screed.openscreed.Open(..., parse_description=True)` 打开文件，根据第一条记录判断格式，选用 FASTA 或 FASTQ 字段类型，再用 `screed.createscreed.create_db` 构建 SQLite 数据库。若输出名不以 screed 的 `_screed` 后缀结尾，则自动补上。随后用 SQL 从数据库读回记录数。
- **API：** `create_screed_database(sequence_file: 'str', output_database: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `screed`（额外依赖 `dnakit[bio-screed]`）；生成本页的环境已安装它，示例是真实运行的。
- **输入：** `sequence_file` 是已存在的 FASTA 或 FASTQ 文件。`output_database` 必须不存在（加上 `_screed` 后缀的名称也不能已存在），其所在文件夹必须存在。
- **输出：** 字典，含 `status`、`backend`（`screed`）、`input_format`（`FASTA` 或 `FASTQ`）、`database`（带后缀的绝对路径）和 `record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import create_screed_database

with tempfile.TemporaryDirectory() as tmp:
    fasta_result = create_screed_database("reads.fa", str(Path(tmp) / "fa_db"), verbose=False)
    print(fasta_result["input_format"], fasta_result["record_count"], Path(fasta_result["database"]).name)
    fastq_result = create_screed_database("reads.fq", str(Path(tmp) / "fq_db_screed"), verbose=False)
    print(fastq_result["input_format"], fastq_result["record_count"], Path(fastq_result["database"]).name)
```

- **示例输出：**

```text
FASTA 2 fa_db_screed
FASTQ 2 fq_db_screed
```

## 23) `dnakit.fasta_fastq.dump_screed_fasta()`

- **作用：** 把 screed 数据库的每条记录写成 FASTA 文本。
- **处理方法：** 检查数据库有记录且含 `name` 和 `sequence` 字段，然后调用 `screed.conversion.ToFasta(database, output)`。标题行是记录名后接描述；来自 FASTQ 的数据库会丢弃质量值。
- **API：** `dump_screed_fasta(database: 'str', output_fasta: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `screed`（额外依赖 `dnakit[bio-screed]`）；生成本页的环境已安装它，示例是真实运行的。
- **输入：** `database` 是已存在的 screed 数据库文件。`output_fasta` 必须不存在，其所在文件夹必须存在。空数据库会抛出 `ConfigurationError`。
- **输出：** 字典，含 `status`、`backend`、`database`、`format`（`FASTA`）、`output`（绝对路径）和 `record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.fasta_fastq import create_screed_database, dump_screed_fasta

with tempfile.TemporaryDirectory() as tmp:
    database = create_screed_database("reads.fq", str(Path(tmp) / "reads_db"), verbose=False)["database"]
    out = Path(tmp) / "reads_dump.fa"
    result = dump_screed_fasta(database, str(out), verbose=False)
    print(result["format"], result["record_count"])
    print(out.read_text(), end="")
```

- **示例输出：**

```text
FASTA 2
>r1 first read
ACGTACGT
>r2 second read
GGCCTTAA
```

## 24) `dnakit.fasta_fastq.dump_screed_fastq()`

- **作用：** 把 screed 数据库的每条记录写成 FASTQ 文本。
- **处理方法：** 检查数据库有记录且含 `name`、`sequence` 和 `quality` 字段，然后调用 `screed.conversion.ToFastq(database, output)`。由 FASTA 构建的数据库没有质量字段，会被 `ConfigurationError` 拒绝。
- **API：** `dump_screed_fastq(database: 'str', output_fastq: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 `screed`（额外依赖 `dnakit[bio-screed]`）；生成本页的环境已安装它，示例是真实运行的。
- **输入：** `database` 是已存在的、由 FASTQ 构建的 screed 数据库。`output_fastq` 必须不存在，其所在文件夹必须存在。
- **输出：** 字典，含 `status`、`backend`、`database`、`format`（`FASTQ`）、`output`（绝对路径）和 `record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.fasta_fastq import create_screed_database, dump_screed_fastq

with tempfile.TemporaryDirectory() as tmp:
    fq_db = create_screed_database("reads.fq", str(Path(tmp) / "fq_db"), verbose=False)["database"]
    out = Path(tmp) / "reads_dump.fq"
    result = dump_screed_fastq(fq_db, str(out), verbose=False)
    print(result["format"], result["record_count"])
    print(out.read_text(), end="")
    fa_db = create_screed_database("reads.fa", str(Path(tmp) / "fa_db"), verbose=False)["database"]
    try:
        dump_screed_fastq(fa_db, str(Path(tmp) / "none.fq"), verbose=False)
    except ConfigurationError as error:
        print(type(error).__name__, "no quality field" in str(error))
```

- **示例输出：**

```text
FASTQ 2
@r1 first read
ACGTACGT
+
IIIIHHHH
@r2 second read
GGCCTTAA
+
########
ConfigurationError True
```
