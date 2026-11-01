# genomepy 与 gget：基因组目录、注释与数据库查询

这一组函数位于 `dnakit.databases`，是对两个可选包和若干外部程序的适配器：[genomepy](https://github.com/vanheeringen-lab/genomepy)（GENCODE、Ensembl、UCSC、NCBI 各 provider 的基因组目录、注释下载与转换、NCBI assembly report、比对软件索引命令）和 [gget](https://github.com/pachterlab/gget)（Ensembl 参考文件链接、NCBI Virus 与 GenBank 下载、cBioPortal 数据与热图、CELLxGENE Census 查询、COSMIC 参考文件，以及 MUSCLE 5 与 DIAMOND 比对程序）。另有五个函数是 genomepy 目录记录的离线变换，完全不需要任何包。gget 函数会在数据旁返回服务名称和访问日期；DNAKit 自身不做缓存（genomepy 会在磁盘上缓存 provider 目录）。

据包元数据，genomepy 使用 MIT 许可，gget 使用 BSD-2-Clause 许可。MUSCLE、DIAMOND、各比对软件和 UCSC 转换工具是各有许可的独立程序；DNAKit 只把它们作为外部进程启动，不复制它们的任何代码。COSMIC 下载函数需要你自己的 COSMIC 账号，账号只用于该次请求。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时函数会抛出带安装命令的 `BackendUnavailableError`；比对软件及 MUSCLE/DIAMOND 的封装在找不到程序时抛出该错误。可用 `dnakit[bio-genomepy]`、`dnakit[bio-gget]` 和 `dnakit[bio-cellxgene]` 安装对应 extra。示例从 `dnakit.databases` 导入。

## 示例数据

本页不声明示例文件。需要输入的示例会把很小的文件写入临时目录；需要互联网、账号或此环境未安装的外部程序的函数，则展示它们对无效输入或程序缺失所抛出的错误。这些小节都附有一个不执行的“安装依赖后的用法”代码块，给出真实调用。

## 1) `dnakit.databases.genomepy_bowtie2_index_command()`

- **作用：** 以参数列表形式运行 `bowtie2-build` 索引命令，并报告新生成的索引文件。
- **处理方法：** 检查 `command` 是列表且第一个元素名为 `bowtie2-build`，从 `executable_path` 或 PATH 解析程序，在 `working_directory` 中不经 shell 运行其余参数（由 DNAKit 的受限运行器约束：单线程、3600 秒超时、1 MB 输出上限），并列出新产生的文件。若程序没有写出任何新文件则抛出 `BackendExecutionError`。它不导入 genomepy 包，而是把 genomepy 的 `bowtie2` 索引插件步骤改写为参数向量形式。
- **API：** `genomepy_bowtie2_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bowtie2-build` 位于 PATH 中，或在 `executable_path` 中给出其路径。不需要 genomepy 包。示例说明：示例传入的命令不以 `bowtie2-build` 开头；检查在查找或运行任何程序之前就失败，因此在任何机器上输出相同。
- **输入：** `command`：非空字符串列表（不接受 shell 字符串），第一个元素的文件名必须是 `bowtie2-build`；`working_directory`：程序运行所在的已存在目录；`executable_path`：程序的可选路径（默认在 PATH 中查找）。
- **输出：** 字典，含 `status`（`complete`）、`operation`（`bowtie2_index_command`）、`backend`（`bowtie2`）、`version`、`outputs`（`index_0`、`index_1` 等，对应 `working_directory` 中新出现的文件路径，递归查找）、`executable`（解析后的程序路径）和 `log_tail`（程序输出的最后 2000 个字符）。
- **示例代码：**

```python
from dnakit.databases import genomepy_bowtie2_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_bowtie2_index_command(
        ["samtools", "faidx", "genome.fa"], working_directory=".", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError The bowtie2 index command must start with 'bowtie2-build'.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_bowtie2_index_command

result = genomepy_bowtie2_index_command(
    ["bowtie2-build", "genome.fa", "genome"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 2) `dnakit.databases.genomepy_bwa_index_command()`

- **作用：** 以参数列表形式运行 `bwa index` 命令，并报告新生成的索引文件。
- **处理方法：** 检查 `command` 是列表且第一个元素名为 `bwa`，从 `executable_path` 或 PATH 解析程序，在 `working_directory` 中不经 shell 运行其余参数（由 DNAKit 的受限运行器约束：单线程、3600 秒超时、1 MB 输出上限），并列出新产生的文件。若程序没有写出任何新文件则抛出 `BackendExecutionError`。它不导入 genomepy 包，而是把 genomepy 的 `bwa` 索引插件步骤改写为参数向量形式。
- **API：** `genomepy_bwa_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `bwa` 位于 PATH 中，或在 `executable_path` 中给出其路径。不需要 genomepy 包。示例说明：此环境未安装 bwa，所以示例写了一个很小的 shell 脚本 `fake_bwa.sh`，它创建两个索引文件，并通过 `executable_path` 传入。函数本身是真实运行的（参数检查、子进程、新文件检测），只是用脚本代替了比对软件。
- **输入：** `command`：非空字符串列表（不接受 shell 字符串），第一个元素的文件名必须是 `bwa`；`working_directory`：程序运行所在的已存在目录；`executable_path`：程序的可选路径（默认在 PATH 中查找）。
- **输出：** 字典，含 `status`（`complete`）、`operation`（`bwa_index_command`）、`backend`（`bwa`）、`version`、`outputs`（`index_0`、`index_1` 等，对应 `working_directory` 中新出现的文件路径，递归查找）、`executable`（解析后的程序路径）和 `log_tail`（程序输出的最后 2000 个字符）。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_bwa_index_command

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("genome.fa").write_text(">chr1\nACGTACGT\n")
    fake = Path("fake_bwa.sh")
    fake.write_text('#!/bin/sh\ntouch "$2.bwt" "$2.sa"\n')
    fake.chmod(0o755)
    result = genomepy_bwa_index_command(
        ["bwa", "index", "genome.fa"],
        working_directory=".",
        executable_path=str(fake.resolve()),
        verbose=False,
    )
    print(result["status"], result["operation"], result["backend"])
    print([Path(path).name for path in result["outputs"].values()])
```

- **示例输出：**

```text
complete bwa_index_command bwa
['genome.fa.bwt', 'genome.fa.sa']
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_bwa_index_command

result = genomepy_bwa_index_command(
    ["bwa", "index", "genome.fa"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 3) `dnakit.databases.genomepy_download_annotation()`

- **作用：** 把注释文件下载（或复制）到基因组文件夹，并像 genomepy 的 base 与 GENCODE provider 那样转换成 GTF 和 BED 文件。
- **处理方法：** 校验参数，检查 UCSC 工具 `genePredToGtf` 在 PATH 中（否则抛出 `BackendUnavailableError`），然后调用 `genomepy.providers.base.download_annotation`。genomepy 会写出 `<genomes_dir>/<localname>/<localname>.annotation.gtf` 和 `.annotation.bed`（转换按输入类型使用 `gtfToGenePred`、`gff3ToGenePred`、`bedToGenePred`、`genePredToGtf`、`genePredToBed`）。随后 DNAKit 列出该文件夹下的全部文件。
- **API：** `genomepy_download_annotation(genomes_dir: 'str', annotation_url: 'str', localname: 'str', *, n: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）、PATH 中的 UCSC 转换工具（例如来自 bioconda），远程 URL 还需要互联网访问。此环境未安装 UCSC 工具，因此示例展示参数校验错误。示例说明：示例使用空的 `localname`，它在检查工具或开始下载之前就被拒绝。
- **输入：** `genomes_dir`：基础文件夹；`annotation_url`：URL（或已存在的本地路径，genomepy 会复制它）；`localname`：非空文件夹名；`n`：可选正整数，只下载开头若干行。文本参数最长 512 个字符。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）。`outputs` 含 `directory` 以及 `file_0`、`file_1` 等（文件夹内全部文件的排序路径）。
- **示例代码：**

```python
from dnakit.databases import genomepy_download_annotation
from dnakit.exceptions import ConfigurationError

try:
    genomepy_download_annotation("genomes", "https://example.org/a.gtf.gz", "", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError localname must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_download_annotation

result = genomepy_download_annotation(
    "genomes",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz",
    "GRCh38_gencode",
)
print(result["outputs"]["directory"])
```

## 4) `dnakit.databases.genomepy_download_annotation_head()`

- **作用：** 只把（可能 gzip 压缩的）注释文件的开头若干行下载到一个新文件。
- **处理方法：** 校验 URL 和新的输出路径（会创建其父文件夹），然后调用 `genomepy.providers.base.download_head`。该函数用 `urllib` 打开 URL，URL 以 `.gz` 结尾时解压，写出每一行，并在写出 `n + 2` 个不以 `#` 开头的行后停止；注释行会写出但不计数。返回的 `lines` 是请求的 `n`，不是实际写出的行数。
- **API：** `genomepy_download_annotation_head(annotation_url: 'str', annotation_file: 'str', *, n: 'int' = 5, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）。远程 URL 需要互联网访问；`file://` URL 可离线使用，示例即采用这种方式。示例说明：示例通过 `file://` URL 读取本地 GTF，因此不访问网络。`n=2` 时文件得到 1 行注释和 4 行特征行（`n + 2`），共 5 行。
- **输入：** `annotation_url`：URL 字符串，最长 512 个字符；`annotation_file`：必须尚不存在的输出路径；`n`：正整数（默认 5）。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）。`outputs` 含 `head`（新文件）；`lines` 是请求的 `n`。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_download_annotation_head

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    rows = ["#!genome-build GRCh38"]
    for i in range(1, 8):
        rows.append(f'chr1\tsrc\tgene\t{i}\t{i + 9}\t.\t+\t.\tgene_id "g{i}";')
    Path("annotation.gtf").write_text("\n".join(rows) + "\n")
    url = Path("annotation.gtf").resolve().as_uri()
    result = genomepy_download_annotation_head(url, "head/head.gtf", n=2, verbose=False)
    print(result["operation"], result["lines"], list(result["outputs"]))
    print(len(Path("head/head.gtf").read_text().splitlines()))
```

- **示例输出：**

```text
annotation_head 2 ['head']
5
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_download_annotation_head

result = genomepy_download_annotation_head(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz",
    "head/gencode_head.gtf",
    n=5,
)
print(result["outputs"]["head"])
```

## 5) `dnakit.databases.genomepy_download_assembly_report()`

- **作用：** 获取某个组装登录号的 NCBI assembly report（染色体名称、角色和长度）。
- **处理方法：** 调用 `genomepy.providers.ncbi.download_assembly_report`。不给 `output_path` 时返回报告各行；给出 `output_path`（新路径）时 genomepy 把报告写成 TSV 文件。登录号不以 `GCA` 或 `GCF` 开头、或 NCBI 没有该组装时，上游不发任何请求就返回 `None`，DNAKit 将其转成 `BackendExecutionError`。DNAKit 自身只检查登录号是非空字符串。
- **API：** `genomepy_download_assembly_report(accession: 'str', *, output_path: 'str | None' = None, quiet: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和互联网访问。示例使用上游在发出任何请求之前就拒绝的登录号，因此可离线运行。
- **输入：** `accession`：字符串，通常为 `GCA_...` 或 `GCF_...`；`output_path`：可选的新文件路径；`quiet`：是否静默 genomepy 警告（默认 True）。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`accession`，以及 `rows`（含报告各列的字典列表）和 `row_count`；若给出 `output_path`，则 `outputs` 含 `report`。
- **示例代码：**

```python
from dnakit.databases import genomepy_download_assembly_report
from dnakit.exceptions import BackendExecutionError

try:
    genomepy_download_assembly_report("XYZ_1", verbose=False)
except BackendExecutionError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
BackendExecutionError No NCBI assembly report for XYZ_1.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_download_assembly_report

result = genomepy_download_assembly_report("GCF_000001405.40")
print(result["row_count"], result["rows"][0]["Sequence-Name"])
```

## 6) `dnakit.databases.genomepy_ensembl_genomes()`

- **作用：** 通过 Ensembl REST API 构建 Ensembl 基因组目录。
- **处理方法：** 调用 `genomepy.providers.ensembl.get_genomes(rest_url)`，它向 REST 服务请求各个 division 及其基因组（跳过 EnsemblBacteria），每个基因组只保留少数摘要字段，并加入固定的 GRCh37 记录。genomepy 会把结果缓存在磁盘上。
- **API：** `genomepy_ensembl_genomes(rest_url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和互联网访问。示例说明：示例传入空 URL，它在发出任何请求之前就被拒绝。
- **输入：** `rest_url`：非空 URL 字符串，最长 512 个字符，例如 `https://rest.ensembl.org/`。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`genomes`（组装名到记录的映射）和 `genome_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_ensembl_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_ensembl_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError rest_url must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_ensembl_genomes

result = genomepy_ensembl_genomes("https://rest.ensembl.org/")
print(result["genome_count"], "GRCh37" in result["genomes"])
```

## 7) `dnakit.databases.genomepy_find_annotation_urls()`

- **作用：** 在存放基因组文件的远程文件夹中查找 GTF 或 GFF3 注释文件。
- **处理方法：** 调用 `genomepy.providers.url.search_url_for_annotations(url)`：取基因组 URL 的文件夹部分，读取该文件夹列表，并按基因组本地名对文件名做模糊匹配。返回匹配文件的链接（没有匹配时为空列表）。
- **API：** `genomepy_find_annotation_urls(url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和互联网访问。示例说明：示例传入空 URL，它在发出任何请求之前就被拒绝。
- **输入：** `url`：基因组文件的非空 URL 字符串（最长 512 个字符）。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`urls`（注释链接列表）和 `url_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_find_annotation_urls
from dnakit.exceptions import ConfigurationError

try:
    genomepy_find_annotation_urls("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError url must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_find_annotation_urls

result = genomepy_find_annotation_urls(
    "https://ftp.ensembl.org/pub/release-112/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz"
)
print(result["urls"])
```

## 8) `dnakit.databases.genomepy_gencode_genomes()`

- **作用：** 构建 GENCODE 基因组目录：每个主要组装的最新发布版本。
- **处理方法：** 调用 `genomepy.providers.gencode.get_genomes(ftp_link)`。它连接 FTP 站点，列出人和小鼠的发布文件夹，按从高到低的顺序遍历发布版本（与 `genomepy_gencode_releases` 的发布规则相同），在每个发布中找到 primary assembly 文件，并为每个组装记录首个找到的发布版本，最后为人类加入 GRCh37 记录。
- **API：** `genomepy_gencode_genomes(ftp_link: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）以及对 GENCODE 服务器的 FTP 访问。示例说明：示例传入空链接，它在建立任何连接之前就被拒绝。
- **输入：** `ftp_link`：非空 FTP 链接字符串，最长 512 个字符，例如 `ftp://ftp.ebi.ac.uk/pub/databases/gencode`。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`genomes`（组装名到含 `annotations`、`taxonomy_id`、`species`、`text_search` 的记录的映射）和 `genome_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_gencode_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_gencode_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError ftp_link must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_gencode_genomes

result = genomepy_gencode_genomes("ftp://ftp.ebi.ac.uk/pub/databases/gencode")
print(sorted(result["genomes"]))
```

## 9) `dnakit.databases.genomepy_gencode_releases()`

- **作用：** 从发布文件夹列表中提取 GENCODE 发布号，按从高到低排序。
- **处理方法：** 检查参数并调用 `genomepy.providers.gencode.get_releases(listing, species)`。对每个末两位字符是数字的条目取其数字；小于等于 21 的丢弃，其余从高到低排序，`mouse` 时在每个前面加 `M`。`README.txt`、`release_M9` 这类条目会被忽略。该函数是纯函数，不发出任何请求。
- **API：** `genomepy_gencode_releases(listing: 'Sequence[str]', species: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 调用的辅助函数来自可选包 genomepy（extra `dnakit[bio-genomepy]`）；不需要互联网访问。
- **输入：** `listing`：字符串列表（单个字符串会被拒绝），通常是以 `release_NN` 或 `release_MNN` 结尾的路径；`species`：`human` 或 `mouse`。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`releases`（字符串列表，从高到低）和 `release_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_gencode_releases

human = genomepy_gencode_releases(
    ["release_44", "release_45", "release_19", "README.txt", "release_43"],
    "human",
    verbose=False,
)
print(human["releases"], human["release_count"])
mouse = genomepy_gencode_releases(
    ["release_M33", "release_M34", "release_M9"], "mouse", verbose=False
)
print(mouse["releases"])
```

- **示例输出：**

```text
['45', '44', '43'] 3
['M34', 'M33']
```

## 10) `dnakit.databases.genomepy_gencode_to_ucsc_names()`

- **作用：** 为一个目录把 GENCODE（Ensembl 风格）的组装名映射成 UCSC 名称。
- **处理方法：** 调用 `genomepy.providers.gencode.get_gencode2ucsc(genomes)`。映射以 `GRCh37` 到 `hg19`、`GRCm38` 到 `mm10` 开始；其他组装按 species 字段决定前缀（以 `H` 开头为 `hg`，否则为 `mm`），再取组装名中的数字作为编号。不发出任何请求。记录缺少 `species` 字段时 genomepy 会抛出 `KeyError`，DNAKit 不做转换。
- **API：** `genomepy_gencode_to_ucsc_names(genomes: 'Mapping[str, Any]', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）；不需要互联网访问。
- **输入：** `genomes`：非空映射，组装名对应带有 `species` 字符串的记录，例如 `genomepy_gencode_genomes` 返回的目录。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`mapping`（GENCODE 名到 UCSC 名）和 `entry_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_gencode_to_ucsc_names

catalogue = {
    "GRCh38": {"species": "Homo sapiens"},
    "GRCm39": {"species": "Mus musculus"},
    "GRCh37": {"species": "Homo sapiens"},
}
result = genomepy_gencode_to_ucsc_names(catalogue, verbose=False)
print(result["mapping"], result["entry_count"])
```

- **示例输出：**

```text
{'GRCh37': 'hg19', 'GRCm38': 'mm10', 'GRCh38': 'hg38', 'GRCm39': 'mm39'} 4
```

## 11) `dnakit.databases.genomepy_generate_annotation()`

- **作用：** 用 UCSC 工具把 GTF 注释转换成 BED，或把 BED 注释转换成 GTF。
- **处理方法：** 要求 `template` 是已存在的文件，且 `target` 已存在时拒绝执行（除非 `overwrite` 为 True）；会创建 `target` 所在文件夹。检查 UCSC 工具 `genePredToGtf` 在 PATH 中（否则 `BackendUnavailableError`），然后调用 `genomepy.annotation.utils.generate_annot`，它把模板依次交给 `gtfToGenePred` 与 `genePredToBed`（GTF 转 BED）或 `bedToGenePred` 与 `genePredToGtf`（BED 转 GTF）。文件名中既无 `gtf` 也无 `bed` 的模板报 `ConfigurationError`；输出文件缺失报 `BackendExecutionError`。
- **API：** `genomepy_generate_annotation(template: 'str', target: 'str', *, overwrite: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和 PATH 中的 UCSC 转换工具（例如 `conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed ucsc-genepredtogtf ucsc-bedtogenepred`）。不需要互联网访问。此环境未安装 UCSC 工具，所以示例展示输入检查。示例说明：示例给出不存在的模板；消息中冒号后面是与机器相关的路径，所以只打印冒号之前的部分。
- **输入：** `template`：已存在的 GTF 或 BED 文件路径；`target`：输出路径；`overwrite`：允许覆盖 `target`（默认 False）。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）。`outputs` 含 `annotation`（转换后的文件）。
- **示例代码：**

```python
from dnakit.databases import genomepy_generate_annotation
from dnakit.exceptions import ConfigurationError

try:
    genomepy_generate_annotation("no_such_file.gtf", "out.bed", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1].split(":")[0])
```

- **示例输出：**

```text
ConfigurationError template does not exist
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_generate_annotation

result = genomepy_generate_annotation("annotation.gtf", "annotation.bed")
print(result["outputs"]["annotation"])
```

## 12) `dnakit.databases.genomepy_genome_size()`

- **作用：** 返回某个 NCBI 组装的总序列长度，由其 assembly report 求和得到。
- **处理方法：** 调用 `genomepy.providers.ncbi.get_genome_size(accession)`，它下载 assembly report 并对 `Sequence-Length` 列求和；ncbi、gencode、ucsc 三个 provider 共用这一实现。没有报告时（例如登录号不以 `GCA` 或 `GCF` 开头，这在不发任何请求的情况下就能判定）上游返回 -1，DNAKit 将其报告为 `genome_size` 为 `None`。
- **API：** `genomepy_genome_size(accession: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和互联网访问。示例使用上游在发出任何请求之前就拒绝的登录号，因此可离线运行。
- **输入：** `accession`：非空字符串，最长 512 个字符，通常为 `GCA_...` 或 `GCF_...`。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`accession` 和 `genome_size`（整数，没有报告时为 `None`）。
- **示例代码：**

```python
from dnakit.databases import genomepy_genome_size

result = genomepy_genome_size("XYZ_1", verbose=False)
print(result["operation"], result["accession"], result["genome_size"])
```

- **示例输出：**

```text
genome_size XYZ_1 None
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_genome_size

result = genomepy_genome_size("GCF_000001405.40")
print(result["genome_size"])
```

## 13) `dnakit.databases.genomepy_gmap_index_command()`

- **作用：** 以参数列表形式运行 `gmap_build` 索引命令，并报告新生成的索引文件。
- **处理方法：** 检查 `command` 是列表且第一个元素名为 `gmap_build`，从 `executable_path` 或 PATH 解析程序，在 `working_directory` 中不经 shell 运行其余参数（由 DNAKit 的受限运行器约束：单线程、3600 秒超时、1 MB 输出上限），并列出新产生的文件。若程序没有写出任何新文件则抛出 `BackendExecutionError`。它不导入 genomepy 包，而是把 genomepy 的 `gmap` 索引插件步骤改写为参数向量形式。
- **API：** `genomepy_gmap_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `gmap_build` 位于 PATH 中，或在 `executable_path` 中给出其路径。不需要 genomepy 包。示例说明：示例把 `working_directory` 指向不存在的文件夹；检查在查找程序之前就失败。冒号之后是与机器相关的路径，所以被截掉。
- **输入：** `command`：非空字符串列表（不接受 shell 字符串），第一个元素的文件名必须是 `gmap_build`；`working_directory`：程序运行所在的已存在目录；`executable_path`：程序的可选路径（默认在 PATH 中查找）。
- **输出：** 字典，含 `status`（`complete`）、`operation`（`gmap_index_command`）、`backend`（`gmap`）、`version`、`outputs`（`index_0`、`index_1` 等，对应 `working_directory` 中新出现的文件路径，递归查找）、`executable`（解析后的程序路径）和 `log_tail`（程序输出的最后 2000 个字符）。
- **示例代码：**

```python
from dnakit.databases import genomepy_gmap_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_gmap_index_command(
        ["gmap_build", "-d", "genome", "genome.fa"],
        working_directory="no_such_directory",
        verbose=False,
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1].split(":")[0])
```

- **示例输出：**

```text
ConfigurationError working_directory does not exist
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_gmap_index_command

result = genomepy_gmap_index_command(
    ["gmap_build", "-D", "indexes", "-d", "genome", "genome.fa"],
    working_directory="indexes",
)
print(sorted(result["outputs"]))
```

## 14) `dnakit.databases.genomepy_head_annotations()`

- **作用：** 对提供方提供的每种注释类型，显示某个基因组注释的前几行。
- **处理方法：** 校验 `name` 与 `lines`，检查 UCSC 工具 `genePredToGtf` 在 PATH 中（否则 `BackendUnavailableError`），并在捕获标准输出的情况下运行 `genomepy.functions.head_annotations(name, provider, n=lines)`。genomepy 遍历列出该基因组的在线 provider，把每种注释的前 `lines` 个基因下载到临时文件夹、转换成 GTF，并打印最多 `lines` 个非注释的 GTF 行。返回其中非空的输出行。
- **API：** `genomepy_head_annotations(name: 'str', *, provider: 'str | None' = None, lines: 'int' = 2, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）、PATH 中的 UCSC 转换工具以及互联网访问。此环境未安装 UCSC 工具，所以示例展示输入检查。示例说明：示例传入空的基因组名，它在工具检查或下载之前就被拒绝。
- **输入：** `name`：基因组名，例如 `hg38`；`provider`：可选的 provider 名；`lines`：正整数（默认 2）。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`genome` 和 `lines`（打印出的行列表）。
- **示例代码：**

```python
from dnakit.databases import genomepy_head_annotations
from dnakit.exceptions import ConfigurationError

try:
    genomepy_head_annotations("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError name must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_head_annotations

result = genomepy_head_annotations("hg38", provider="UCSC", lines=2)
print(result["lines"])
```

## 15) `dnakit.databases.genomepy_hisat2_index_command()`

- **作用：** 以参数列表形式运行 `hisat2-build` 索引命令，并报告新生成的索引文件。
- **处理方法：** 检查 `command` 是列表且第一个元素名为 `hisat2-build`，从 `executable_path` 或 PATH 解析程序，在 `working_directory` 中不经 shell 运行其余参数（由 DNAKit 的受限运行器约束：单线程、3600 秒超时、1 MB 输出上限），并列出新产生的文件。若程序没有写出任何新文件则抛出 `BackendExecutionError`。它不导入 genomepy 包，而是把 genomepy 的 `hisat2` 索引插件步骤改写为参数向量形式。
- **API：** `genomepy_hisat2_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `hisat2-build` 位于 PATH 中，或在 `executable_path` 中给出其路径。不需要 genomepy 包。示例说明：示例把命令作为单个 shell 字符串传入，函数会拒绝；不会查找或运行任何程序。
- **输入：** `command`：非空字符串列表（不接受 shell 字符串），第一个元素的文件名必须是 `hisat2-build`；`working_directory`：程序运行所在的已存在目录；`executable_path`：程序的可选路径（默认在 PATH 中查找）。
- **输出：** 字典，含 `status`（`complete`）、`operation`（`hisat2_index_command`）、`backend`（`hisat2`）、`version`、`outputs`（`index_0`、`index_1` 等，对应 `working_directory` 中新出现的文件路径，递归查找）、`executable`（解析后的程序路径）和 `log_tail`（程序输出的最后 2000 个字符）。
- **示例代码：**

```python
from dnakit.databases import genomepy_hisat2_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_hisat2_index_command(
        "hisat2-build genome.fa genome", working_directory=".", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError command must be a non-empty list of strings (no shell string).
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_hisat2_index_command

result = genomepy_hisat2_index_command(
    ["hisat2-build", "genome.fa", "genome"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 16) `dnakit.databases.genomepy_list_available_genomes()`

- **作用：** 列出 genomepy 在线 provider 提供的基因组。
- **处理方法：** 调用 `genomepy.functions.list_available_genomes(provider, size)`，向各个在线 provider（GENCODE、Ensembl、UCSC、NCBI）请求其目录，并把各行转换成普通列表。每行为 `[name, provider, description...]`。DNAKit 不做额外校验：provider 名未知时 genomepy 抛出 `ValueError`，不做转换。
- **API：** `genomepy_list_available_genomes(*, provider: 'str | None' = None, include_size: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和互联网访问。示例使用未知的 provider 名，genomepy 在发出任何请求之前就会拒绝。
- **输入：** `provider`：可选的 provider 名（默认全部在线 provider）；`include_size`：在每行中加入基因组绝对大小（默认 False）。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`rows`（列表的列表）和 `row_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_list_available_genomes

try:
    genomepy_list_available_genomes(provider="nonexistent_provider", verbose=False)
except ValueError as error:
    print(type(error).__name__, str(error).split(".")[0])
```

- **示例输出：**

```text
ValueError Unknown provider 'nonexistent_provider'
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_list_available_genomes

result = genomepy_list_available_genomes(provider="UCSC")
print(result["row_count"], result["rows"][0][:2])
```

## 17) `dnakit.databases.genomepy_minimap2_index_command()`

- **作用：** 以参数列表形式运行 `minimap2` 索引命令，并报告新生成的索引文件。
- **处理方法：** 检查 `command` 是列表且第一个元素名为 `minimap2`，从 `executable_path` 或 PATH 解析程序，在 `working_directory` 中不经 shell 运行其余参数（由 DNAKit 的受限运行器约束：单线程、3600 秒超时、1 MB 输出上限），并列出新产生的文件。若程序没有写出任何新文件则抛出 `BackendExecutionError`。它不导入 genomepy 包，而是把 genomepy 的 `minimap2` 索引插件步骤改写为参数向量形式。
- **API：** `genomepy_minimap2_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `minimap2` 位于 PATH 中，或在 `executable_path` 中给出其路径。不需要 genomepy 包。示例说明：示例传入空的命令列表，函数会拒绝；不会查找或运行任何程序。
- **输入：** `command`：非空字符串列表（不接受 shell 字符串），第一个元素的文件名必须是 `minimap2`；`working_directory`：程序运行所在的已存在目录；`executable_path`：程序的可选路径（默认在 PATH 中查找）。
- **输出：** 字典，含 `status`（`complete`）、`operation`（`minimap2_index_command`）、`backend`（`minimap2`）、`version`、`outputs`（`index_0`、`index_1` 等，对应 `working_directory` 中新出现的文件路径，递归查找）、`executable`（解析后的程序路径）和 `log_tail`（程序输出的最后 2000 个字符）。
- **示例代码：**

```python
from dnakit.databases import genomepy_minimap2_index_command
from dnakit.exceptions import ConfigurationError

try:
    genomepy_minimap2_index_command([], working_directory=".", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError command must be a non-empty list of strings (no shell string).
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_minimap2_index_command

result = genomepy_minimap2_index_command(
    ["minimap2", "-d", "genome.mmi", "genome.fa"], working_directory="indexes"
)
print(sorted(result["outputs"]))
```

## 18) `dnakit.databases.genomepy_ncbi_genomes()`

- **作用：** 把 NCBI assembly summary 文件解析成基因组目录。
- **处理方法：** 调用 `genomepy.providers.ncbi.get_genomes(assembly_url)`。genomepy 用 `urllib` 打开该 URL 下的四个文件（`assembly_summary_genbank_historical.txt`、`assembly_summary_refseq_historical.txt`、`assembly_summary_genbank.txt`、`assembly_summary_refseq.txt`），跳过注释行、读取表头行，并保留每行的八列。`ftp_path` 为 `na` 的行被丢弃。后面的文件会覆盖前面同名的组装（第 15 列，经安全处理后作为名称）。
- **API：** `genomepy_ncbi_genomes(assembly_url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）。真实的 NCBI 摘要文件很大，需要互联网访问；指向本地文件的 `file://` URL 可离线使用，示例即采用这种方式。示例说明：示例写入四个很小的摘要文件（RefSeq 文件中含一个组装），并把它们所在文件夹作为 `file://` URL 传入，因此不访问网络。
- **输入：** `assembly_url`：存放四个摘要文件的文件夹的非空 URL 字符串（最长 512 个字符），例如 `https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS`。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`genomes`（组装名到含 `assembly_accession`、`taxid`、`species_taxid`、`organism_name`、`asm_submitter`、`gbrs_paired_asm`、`paired_asm_comp`、`ftp_path` 的记录的映射）和 `genome_count`。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_ncbi_genomes

COLUMNS = [
    "assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category",
    "taxid", "species_taxid", "organism_name", "infraspecific_name", "isolate",
    "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date",
    "asm_name", "asm_submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
]
FILES = [
    "assembly_summary_genbank_historical.txt",
    "assembly_summary_refseq_historical.txt",
    "assembly_summary_genbank.txt",
    "assembly_summary_refseq.txt",
]

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("summaries").mkdir()
    record = dict.fromkeys(COLUMNS, "na")
    record.update(
        assembly_accession="GCF_000001405.40", taxid="9606", species_taxid="9606",
        organism_name="Homo sapiens", asm_name="GRCh38.p14", asm_submitter="Genome Reference Consortium",
        gbrs_paired_asm="GCA_000001405.29", paired_asm_comp="identical",
        ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
    )
    for name in FILES:
        text = "# comment line\n#" + "\t".join(COLUMNS) + "\n"
        if name == "assembly_summary_refseq.txt":
            text += "\t".join(record[column] for column in COLUMNS) + "\n"
        Path("summaries", name).write_text(text)
    url = Path("summaries").resolve().as_uri()
    result = genomepy_ncbi_genomes(url, verbose=False)
    genome = result["genomes"]["GRCh38.p14"]
    print(result["genome_count"], list(result["genomes"]))
    print(genome["assembly_accession"], genome["organism_name"], genome["taxid"])
```

- **示例输出：**

```text
1 ['GRCh38.p14']
GCF_000001405.40 Homo sapiens 9606
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_ncbi_genomes

result = genomepy_ncbi_genomes("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS")
print(result["genome_count"])
```

## 19) `dnakit.databases.genomepy_run_index_cmd()`

- **作用：** 按名称选择一个比对软件（bowtie2、bwa、gmap、hisat2、minimap2 或 star），以参数列表形式运行其索引命令。
- **处理方法：** 检查 `name` 属于六个比对软件之一，然后执行与对应的 `genomepy_<aligner>_index_command` 封装相同的步骤：`command` 的第一个元素必须以该软件的索引程序命名（`bowtie2-build`、`bwa`、`gmap_build`、`hisat2-build`、`minimap2`、`STAR`），程序从 `executable_path` 或 PATH 解析，在 `working_directory` 中不经 shell 运行，并列出新文件。与把字符串交给 shell 的 genomepy 原生 `run_index_cmd` 不同，它只接受参数列表。
- **API：** `genomepy_run_index_cmd(name: 'str', command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要所选比对软件的外部程序位于 PATH 中，或在 `executable_path` 中给出其路径。不需要 genomepy 包。示例说明：示例给出不在列表中的比对软件名；检查在查找任何程序之前就失败。
- **输入：** `name`：`bowtie2`、`bwa`、`gmap`、`hisat2`、`minimap2`、`star` 之一；`command`、`working_directory`、`executable_path` 与专用封装相同。
- **输出：** 字典，含 `status`（`complete`）、`operation`（`<name>_index_command`）、`backend`（比对软件名）、`version`、`outputs`（新文件对应的 `index_0` 等）、`executable` 和 `log_tail`。
- **示例代码：**

```python
from dnakit.databases import genomepy_run_index_cmd
from dnakit.exceptions import ConfigurationError

try:
    genomepy_run_index_cmd("blast", ["makeblastdb"], working_directory=".", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError name must be one of ['bowtie2', 'bwa', 'gmap', 'hisat2', 'minimap2', 'star'].
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_run_index_cmd

result = genomepy_run_index_cmd(
    "star",
    ["STAR", "--runMode", "genomeGenerate", "--genomeDir", "star_index",
     "--genomeFastaFiles", "genome.fa"],
    working_directory="indexes",
)
print(sorted(result["outputs"]))
```

## 20) `dnakit.databases.genomepy_search_genomes()`

- **作用：** 在各在线 provider 中按基因组名称、描述、登录号和分类 ID 搜索。
- **处理方法：** 调用 `genomepy.providers.search(term, provider, exact, size)`。genomepy 对名称和描述做不区分大小写的搜索（允许正则表达式；多个词不分顺序），以 `GCA_` 或 `GCF_` 开头的词与组装登录号匹配，数字与分类 ID 匹配。每行为 `[name, provider, description...]`。
- **API：** `genomepy_search_genomes(term: 'str', *, provider: 'str | None' = None, exact: 'bool' = False, include_size: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和互联网访问。示例说明：示例传入空的搜索词，它在联系任何 provider 之前就被拒绝。
- **输入：** `term`：非空字符串，最长 512 个字符；`provider`：可选的 provider 名；`exact`：要求精确匹配（默认 False）；`include_size`：在每行中加入基因组大小（默认 False）。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`rows`（列表的列表）和 `row_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_search_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_search_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError term must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_search_genomes

result = genomepy_search_genomes("zebrafish", provider="Ensembl")
print(result["row_count"], result["rows"][0][:2])
```

## 21) `dnakit.databases.genomepy_star_index_command()`

- **作用：** 以参数列表形式运行 STAR 的基因组索引生成命令，并报告新生成的索引文件。
- **处理方法：** 检查 `command` 是列表且第一个元素名为 `STAR`，从 `executable_path` 或 PATH 解析程序，在 `working_directory` 中不经 shell 运行其余参数（由 DNAKit 的受限运行器约束：单线程、3600 秒超时、1 MB 输出上限），并列出新产生的文件。若程序没有写出任何新文件则抛出 `BackendExecutionError`。它不导入 genomepy 包，而是把 genomepy 的 `star` 索引插件步骤改写为参数向量形式。
- **API：** `genomepy_star_index_command(command: 'Sequence[str]', *, working_directory: 'str', executable_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 `STAR` 位于 PATH 中，或在 `executable_path` 中给出其路径。不需要 genomepy 包。示例说明：此环境未安装 STAR，所以示例写了一个很小的 shell 脚本 `fake_star.sh`，它创建三个文件，并通过 `executable_path` 传入。函数本身是真实运行的，只是用脚本代替了比对软件。
- **输入：** `command`：非空字符串列表（不接受 shell 字符串），第一个元素的文件名必须是 `STAR`；`working_directory`：程序运行所在的已存在目录；`executable_path`：程序的可选路径（默认在 PATH 中查找）。
- **输出：** 字典，含 `status`（`complete`）、`operation`（`star_index_command`）、`backend`（`star`）、`version`、`outputs`（`index_0`、`index_1` 等，对应 `working_directory` 中新出现的文件路径，递归查找）、`executable`（解析后的程序路径）和 `log_tail`（程序输出的最后 2000 个字符）。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import genomepy_star_index_command

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("genome.fa").write_text(">chr1\nACGTACGT\n")
    fake = Path("fake_star.sh")
    fake.write_text('#!/bin/sh\ntouch Genome SA SAindex\n')
    fake.chmod(0o755)
    result = genomepy_star_index_command(
        ["STAR", "--runMode", "genomeGenerate", "--genomeFastaFiles", "genome.fa"],
        working_directory=".",
        executable_path=str(fake.resolve()),
        verbose=False,
    )
    print(result["operation"], result["backend"])
    print([Path(path).name for path in result["outputs"].values()])
```

- **示例输出：**

```text
star_index_command star
['Genome', 'SA', 'SAindex']
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_star_index_command

result = genomepy_star_index_command(
    ["STAR", "--runMode", "genomeGenerate", "--genomeDir", "star_index",
     "--genomeFastaFiles", "genome.fa"],
    working_directory="indexes",
)
print(sorted(result["outputs"]))
```

## 22) `dnakit.databases.genomepy_ucsc_download_annotation()`

- **作用：** 从 UCSC MySQL 服务器下载 UCSC 基因注释（扩展 genePred），并转换成 GTF 和 BED 文件。
- **处理方法：** 校验文本参数，检查 UCSC 工具 `genePredToGtf` 在 PATH 中（否则 `BackendUnavailableError`），然后调用 `genomepy.providers.ucsc.download_annotation`。genomepy 读取注释表的列，从 MySQL 服务器查询表内容，写出扩展 genePred 文件，并用 `genePredToGtf` 和 `genePredToBed` 转换到 `<genomes_dir>/<localname>/`。DNAKit 列出该文件夹下的所有文件。
- **API：** `genomepy_ucsc_download_annotation(name: 'str', annotation: 'str', genomes_dir: 'str', localname: 'str', *, n: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）、PATH 中的 UCSC 转换工具，以及对 UCSC MySQL 服务器的互联网访问。此环境未安装 UCSC 工具，所以示例展示输入检查。示例说明：示例传入空的基因组名，它在工具检查或下载之前就被拒绝。
- **输入：** `name`：UCSC 基因组（数据库）名，例如 `hg38`；`annotation`：表名，例如 `ncbiRefSeq`、`refGene`、`ensGene` 或 `knownGene`；`genomes_dir`：基础文件夹；`localname`：输出文件夹名；`n`：可选正整数，限制行数。文本参数非空且最长 512 个字符。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）。`outputs` 含 `directory` 以及 `file_0`、`file_1` 等（文件夹内文件的排序路径）。
- **示例代码：**

```python
from dnakit.databases import genomepy_ucsc_download_annotation
from dnakit.exceptions import ConfigurationError

try:
    genomepy_ucsc_download_annotation("", "refGene", "genomes", "hg38", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError name must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_ucsc_download_annotation

result = genomepy_ucsc_download_annotation("hg38", "refGene", "genomes", "hg38_refgene")
print(result["outputs"]["directory"])
```

## 23) `dnakit.databases.genomepy_ucsc_genomes()`

- **作用：** 通过 UCSC REST API 构建 UCSC 基因组目录。
- **处理方法：** 调用 `genomepy.providers.ucsc.get_genomes(rest_url)`，它请求 UCSC 基因组的 JSON 列表，每个基因组保留五个字段（`taxId`、`description`、`scientificName`、`sourceName`、`htmlPath`），并加入在 UCSC MySQL 服务器上找到的注释表和组装登录号（这些记录变换以 `genomepy_ucsc_assembly_accessions`、`genomepy_ucsc_sourcename_accessions` 和 `genomepy_ucsc_annotation_links` 的形式提供离线版本）。
- **API：** `genomepy_ucsc_genomes(rest_url: 'str', *, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 genomepy（extra `dnakit[bio-genomepy]`）和互联网访问。示例说明：示例传入空 URL，它在发出任何请求之前就被拒绝。
- **输入：** `rest_url`：非空 URL 字符串，最长 512 个字符，例如 `https://api.genome.ucsc.edu/list/ucscGenomes`。
- **输出：** 字典，含genomepy 通用结果字段 `status`（`complete`）、`operation`、`backend`、`version`（已安装的 genomepy 版本，缺失时为 `None`）和 `outputs`（名称到路径字符串的映射）、`genomes`（UCSC 组装名到其记录的映射）和 `genome_count`。
- **示例代码：**

```python
from dnakit.databases import genomepy_ucsc_genomes
from dnakit.exceptions import ConfigurationError

try:
    genomepy_ucsc_genomes("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError rest_url must be a non-empty string of at most 512 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import genomepy_ucsc_genomes

result = genomepy_ucsc_genomes("https://api.genome.ucsc.edu/list/ucscGenomes")
print(result["genome_count"], "hg38" in result["genomes"])
```

## 24) `dnakit.databases.genomepy_ucsc_assembly_accessions()`

- **作用：** 根据 UCSC 的组装等价关系行，为 UCSC 目录设置 `assembly_accession`。
- **处理方法：** 这是 genomepy UCSC provider 中 `add_accessions1` 步骤的离线版本，上游从 UCSC MySQL 服务器读取 `asmEquivalent` 表，这里改为直接传入各行。对每个源组装只保留匹配数最高的行；从目标名中截取 `GCA_`/`GCF_` 登录号；某组装保留多行时只使用 `GCA_` 登录号。不在 `genomes` 中的源会被忽略，保留行的目标名不含登录号时不赋值（上游会存入 NaN）。输入会被复制，不会被修改。
- **API：** `genomepy_ucsc_assembly_accessions(genomes: 'Mapping[str, Mapping[str, Any]]', equivalences: 'Sequence[Sequence[Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`。
- **运行要求：** 不需要可选包，也不需要互联网访问；普通 pip 安装后即可使用。示例说明：`hg19` 有两行并列最高匹配数 900，所以只使用 `GCA_` 登录号；`zzz` 不在目录中，被忽略；输入目录保持不变。
- **输入：** `genomes`：UCSC 组装名到记录（映射）的映射；`equivalences`：由恰好三个字段 `(source, destination, matchCount)` 组成的行。类型错误会抛出 `ConfigurationError`。
- **输出：** 新的目录（字典），键与输入相同；获得登录号的记录带有字符串 `assembly_accession`。
- **示例代码：**

```python
from dnakit.databases import genomepy_ucsc_assembly_accessions

catalogue = {"hg19": {}, "hg38": {}, "mm10": {}}
rows = [
    ("hg19", "GCA_000001405.1", 900),
    ("hg19", "GCF_000001405.13", 900),
    ("hg19", "GCA_000000000.1", 100),
    ("hg38", "GCF_000001405.39", 5000),
    ("mm10", "GCA_000001635.2", 10),
    ("zzz", "GCA_000000001.1", 3),
]
result = genomepy_ucsc_assembly_accessions(catalogue, rows, verbose=False)
for name, record in result.items():
    print(name, record["assembly_accession"])
print(catalogue)
```

- **示例输出：**

```text
hg19 GCA_000001405.1
hg38 GCF_000001405.39
mm10 GCA_000001635.2
{'hg19': {}, 'hg38': {}, 'mm10': {}}
```

## 25) `dnakit.databases.genomepy_ucsc_sourcename_accessions()`

- **作用：** 根据每条记录 `sourceName` 中的登录号，为 UCSC 目录设置 `assembly_accession`。
- **处理方法：** 这是 genomepy UCSC provider 中 `add_accessions2` 步骤的离线版本。对每条记录，保存 `sourceName` 中第一个 `GCA_`/`GCF_` 登录号（九位数字加版本号）。缺少 `sourceName` 的记录抛出 `ConfigurationError`（上游为 `KeyError`）；`sourceName` 为 `None` 时跳过（上游为 `TypeError`）；非字符串值抛出 `ConfigurationError`。输入会被复制，不会被修改。
- **API：** `genomepy_ucsc_sourcename_accessions(genomes: 'Mapping[str, Mapping[str, Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`。
- **运行要求：** 不需要可选包，也不需要互联网访问；普通 pip 安装后即可使用。
- **输入：** `genomes`：组装名到记录的映射，记录须含 `sourceName` 字段（字符串或 `None`）。
- **输出：** 新的目录（字典）；`sourceName` 含登录号的记录会增加 `assembly_accession`。
- **示例代码：**

```python
from dnakit.databases import genomepy_ucsc_sourcename_accessions
from dnakit.exceptions import ConfigurationError

catalogue = {
    "mm39": {"sourceName": "Genome Reference Consortium GCA_000001635.9"},
    "hg38": {"sourceName": "GRCh38 GCF_000001405.40"},
    "unknown": {"sourceName": None},
}
result = genomepy_ucsc_sourcename_accessions(catalogue, verbose=False)
for name, record in result.items():
    print(name, record.get("assembly_accession"))
print("assembly_accession" in catalogue["mm39"])
try:
    genomepy_ucsc_sourcename_accessions({"bad": {}}, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
mm39 GCA_000001635.9
hg38 GCF_000001405.40
unknown None
False
ConfigurationError genome record 'bad' has no sourceName field.
```

## 26) `dnakit.databases.genomepy_ucsc_annotation_links()`

- **作用：** 把 UCSC 为每个基因组列出的注释表名追加到其记录的 `annotations` 列表。
- **处理方法：** 这是 genomepy UCSC provider 中 `add_annotation_links` 步骤的离线版本。上游查询 UCSC MySQL 服务器的 `information_schema`，获取表 `ensGene`、`ncbiRefSeq`、`knownGene` 和 `refGene`；这里直接传入该查询结果的 `(database, table)` 行，并按行顺序把每个表名追加到同名数据库对应基因组的 `annotations` 列表。不在目录中的基因组的行会被忽略。匹配的记录没有 `annotations` 列表时抛出 `ConfigurationError`（上游为 `KeyError`）。输入会被复制，不会被修改。
- **API：** `genomepy_ucsc_annotation_links(genomes: 'Mapping[str, Mapping[str, Any]]', table_rows: 'Sequence[Sequence[Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`。
- **运行要求：** 不需要可选包，也不需要互联网访问；普通 pip 安装后即可使用。
- **输入：** `genomes`：UCSC 数据库名到含 `annotations` 列表的记录的映射；`table_rows`：由恰好两个字段 `(database, table)` 组成的行。
- **输出：** 新的目录（字典），其 `annotations` 列表包含这些表名。
- **示例代码：**

```python
from dnakit.databases import genomepy_ucsc_annotation_links

catalogue = {"hg38": {"annotations": []}, "mm39": {"annotations": ["existing"]}}
rows = [("hg38", "ensGene"), ("hg38", "refGene"), ("mm39", "knownGene"), ("other", "ensGene")]
result = genomepy_ucsc_annotation_links(catalogue, rows, verbose=False)
print(result)
print(catalogue)
```

- **示例输出：**

```text
{'hg38': {'annotations': ['ensGene', 'refGene']}, 'mm39': {'annotations': ['existing', 'knownGene']}}
{'hg38': {'annotations': []}, 'mm39': {'annotations': ['existing']}}
```

## 27) `dnakit.databases.genomepy_ensembl_add_grch37()`

- **作用：** 加入 genomepy 在 Ensembl 发布目录之外单独列出的固定 GRCh37 记录。
- **处理方法：** 这是 genomepy Ensembl provider 中 `add_grch37` 的离线版本。记录是固定的字面量（智人、组装登录号 `GCA_000001405.14`、分类 ID 9606、division 为 `vertebrates`、Ensembl GRCh37 第 87 版 GTF），以键 `GRCh37` 插入；已有的 `GRCh37` 条目会被替换。`genome` 与 `annotation` 的 URL 中带有发布号占位符 `{}`，与上游一致。输入会被复制，不会被修改。
- **API：** `genomepy_ensembl_add_grch37(genomes: 'Mapping[str, Mapping[str, Any]]', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`。
- **运行要求：** 不需要可选包，也不需要互联网访问；普通 pip 安装后即可使用。
- **输入：** `genomes`：Ensembl 目录，即组装名到记录（映射）的映射。
- **输出：** 新的目录（字典），多出键 `GRCh37`；其记录含 `genome`、`annotation`、`assembly_accession`、`taxonomy_id`、`name`、`scientific_name`、`url_name`、`assembly_name`、`division`、`base_count`、`display_name` 和 `genebuild`。
- **示例代码：**

```python
from dnakit.databases import genomepy_ensembl_add_grch37

catalogue = {"GRCh38": {"name": "Homo_sapiens"}}
result = genomepy_ensembl_add_grch37(catalogue, verbose=False)
record = result["GRCh37"]
print(list(result), list(catalogue))
print(record["assembly_accession"], record["taxonomy_id"], record["scientific_name"])
print(record["annotation"])
```

- **示例输出：**

```text
['GRCh38', 'GRCh37'] ['GRCh38']
GCA_000001405.14 9606 Homo sapiens
http://ftp.ensembl.org/pub/grch37/release-{}/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
```

## 28) `dnakit.databases.genomepy_gencode_add_grch37()`

- **作用：** 加入 GENCODE 放在专门映射文件夹中的 GRCh37 注释记录。
- **处理方法：** 这是 genomepy GENCODE provider 中 `add_grch37` 的离线版本。在形如 `GRCh<数字>` 的键中取最大者，从该组装第一条注释链接的 `release_` 段读取发布号（取该段最后两个字符），并构造 GRCh37 lift-over 链接 `<ftp_link>/Gencode_human/release_<N>/GRCh37_mapping/gencode.v<N>lift37.annotation.gtf.gz`。与上游对所有键取最大值不同，这里只考虑 `GRCh` 键；目录中没有可用的 `GRCh` 键、注释链接或 `release_` 段时抛出 `ConfigurationError`。输入会被复制，不会被修改。
- **API：** `genomepy_gencode_add_grch37(genomes: 'Mapping[str, Mapping[str, Any]]', ftp_link: 'str', *, verbose: 'bool' = True) -> 'dict[str, dict[str, Any]]'`。
- **运行要求：** 不需要可选包，也不需要互联网访问；普通 pip 安装后即可使用。
- **输入：** `genomes`：GENCODE 目录，即组装名到含 `annotations` 列表的记录的映射；`ftp_link`：GENCODE 列表的非空基础链接，例如 `ftp://ftp.ebi.ac.uk/pub/databases/gencode`。
- **输出：** 新的目录（字典），多出键 `GRCh37`；其记录含 `annotations`（一个 URL）、`taxonomy_id`（9606）、`species` 和 `text_search`。
- **示例代码：**

```python
from dnakit.databases import genomepy_gencode_add_grch37

base = "ftp://ftp.ebi.ac.uk/pub/databases/gencode"
catalogue = {
    "GRCh38": {"annotations": [base + "/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"]},
    "GRCm39": {"annotations": []},
}
result = genomepy_gencode_add_grch37(catalogue, base, verbose=False)
record = result["GRCh37"]
print(list(result), list(catalogue))
print(record["annotations"][0])
print(record["taxonomy_id"], record["species"], record["text_search"])
```

- **示例输出：**

```text
['GRCh38', 'GRCm39', 'GRCh37'] ['GRCh38', 'GRCm39']
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.annotation.gtf.gz
9606 Homo sapiens human
```

## 29) `dnakit.databases.gget_alphainfluenza_download()`

- **作用：** 把甲型流感（alphainfluenza）序列从 NCBI 缓存数据包下载到一个文件夹。
- **处理方法：** 必要时创建 `outdir`，然后用可选筛选条件（宿主、仅完整基因组、是否注释、登录号以及是否使用登录号）调用 `gget.gget_virus.download_alphainfluenza_optimized`。gget 返回 ZIP 路径、NCBI 能够应用的筛选条件，以及无法应用的条件。
- **API：** `gget_alphainfluenza_download(*, outdir: 'str', host: 'str | None' = None, complete_only: 'bool | None' = None, annotated: 'bool | None' = None, accession: 'str | None' = None, use_accession: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入空的文件夹路径，它在导入 gget 或发出任何请求之前就被拒绝。
- **输入：** `outdir`：输出文件夹路径（缺失时创建）；`host`、`accession`：可选字符串；`complete_only`、`annotated`：可选布尔值（`None` 表示不设置该筛选）；`use_accession`：布尔值（默认 False）。非布尔标志会抛出 `ConfigurationError`。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`outdir`、`zip_path`、`applied_filters` 和 `missing_filters`。
- **示例代码：**

```python
from dnakit.databases import gget_alphainfluenza_download
from dnakit.exceptions import ConfigurationError

try:
    gget_alphainfluenza_download(outdir="", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError outdir must be a folder path.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_alphainfluenza_download

result = gget_alphainfluenza_download(outdir="influenza_a", host="human", complete_only=True)
print(result["zip_path"], result["applied_filters"])
```

## 30) `dnakit.databases.gget_cbio_download()`

- **作用：** 把 cBioPortal 研究数据包（突变、拷贝数、结构变异、临床数据）下载到一个文件夹。
- **处理方法：** 校验研究 ID，必要时创建 `out_dir`，然后调用 `gget.gget_cbio.download_cbioportal_data(study_ids, verbose=False, out_dir, confirm_download)`。gget 报告并非所有文件都已下载时，函数抛出 `BackendExecutionError`。DNAKit 的 `verbose` 参数不会传给 gget。
- **API：** `gget_cbio_download(study_ids: 'Sequence[str]', *, out_dir: 'str', confirm_download: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入含空格的研究 ID，它在创建任何文件夹之前就未通过标识符检查。
- **输入：** `study_ids`：1 至 10000 个标识符的列表（字母、数字和 `._:/+-`，每个最长 64 个字符）；`out_dir`：输出文件夹（缺失时创建）；`confirm_download`：传给 gget（默认 False）。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`study_ids` 和 `out_dir`。
- **示例代码：**

```python
from dnakit.databases import gget_cbio_download
from dnakit.exceptions import ConfigurationError

try:
    gget_cbio_download(["bad study id"], out_dir="cbio", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError study_ids[] has an unsupported format: 'bad study id'.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_cbio_download

result = gget_cbio_download(["msk_impact_2017"], out_dir="cbio_data")
print(result["out_dir"])
```

## 31) `dnakit.databases.gget_cbio_plot()`

- **作用：** 下载 cBioPortal 研究数据，并绘制每个分层中基因变异的热图。
- **处理方法：** 校验研究 ID、最多 200 个基因符号、`stratification`、`variation_type`、不含路径分隔符的 `figure_filename` 以及 `dpi`（30 至 1200）；必要时创建 `figure_dir` 和 `data_dir`；然后调用 `gget.cbio_plot(..., show=False)`。gget 返回假值时函数抛出 `BackendExecutionError`。仅当 PNG 之后确实存在时才返回图片路径。
- **API：** `gget_cbio_plot(study_ids: 'Sequence[str]', genes: 'Sequence[str]', *, figure_dir: 'str', figure_filename: 'str', data_dir: 'str', stratification: 'str' = 'tissue', variation_type: 'str' = 'mutation_occurrences', confirm_download: 'bool' = False, dpi: 'int' = 100, figure_title: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入空的研究列表，它在创建任何文件夹之前就被拒绝。
- **输入：** `study_ids`：研究标识符列表；`genes`：1 至 200 个基因符号的列表；`figure_dir`、`figure_filename`、`data_dir`：显式的输出位置；`stratification`（默认 `tissue`）和 `variation_type`（默认 `mutation_occurrences`）：传给 gget；`confirm_download`、`dpi`（默认 100）、`figure_title`。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`study_ids`、`genes`、`figure`（PNG 路径，或 `None`）和 `data_dir`。
- **示例代码：**

```python
from dnakit.databases import gget_cbio_plot
from dnakit.exceptions import ConfigurationError

try:
    gget_cbio_plot(
        [], ["TP53"], figure_dir="figs", figure_filename="a.png", data_dir="data", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError study_ids must hold between 1 and 10000 identifiers.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_cbio_plot

result = gget_cbio_plot(
    ["msk_impact_2017"],
    ["TP53", "KRAS"],
    figure_dir="figures",
    figure_filename="alterations.png",
    data_dir="cbio_data",
)
print(result["figure"])
```

## 32) `dnakit.databases.gget_diamond_align()`

- **作用：** 用 DIAMOND 把查询序列比对到参考序列，并返回表格形式的命中。
- **处理方法：** 检查查询文件与参考文件、新的输出路径、灵敏度名称和线程数，然后查找 `diamond` 程序（来自 `diamond_path` 或 PATH；找不到时抛出带安装提示的 `BackendUnavailableError`）并读取其版本。在 `output_path` 旁的临时文件夹中先对参考序列运行 `diamond makedb`，再运行 `diamond blastp`（`translated` 为 True 时为 `blastx`），使用 `--outfmt 6` 与固定的 14 列，随后删除该文件夹。输出文件被解析为每个命中一个字典。
- **API：** `gget_diamond_align(query: 'str', reference: 'str', output_path: 'str', *, translated: 'bool' = False, sensitivity: 'str' = 'very-sensitive', threads: 'int' = 1, diamond_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 DIAMOND（2.1 或更高）位于 PATH 中，或在 `diamond_path` 中给出其路径（gget 在 `gget/bins/Linux/diamond` 下附带 Linux 版本）。不会导入 gget 包本身。此环境的 PATH 中没有 DIAMOND，所以示例展示程序缺失错误。示例说明：示例把 `diamond_path` 指向不存在的文件，因此无论是否安装 DIAMOND 都会以 `BackendUnavailableError` 失败。
- **输入：** `query`、`reference`：已存在的 FASTA 文件；`output_path`：新的文件路径，其所在文件夹须存在；`translated`：False 为蛋白对蛋白，True 为核酸对蛋白（默认 False）；`sensitivity`：`fast`、`mid-sensitive`、`sensitive`、`more-sensitive`、`very-sensitive` 或 `ultra-sensitive`（默认 `very-sensitive`）；`threads`：1 至 256；`diamond_path`：可选的程序路径。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`program`（名称、版本、路径）、`mode`（`blastp` 或 `blastx`）、`sensitivity`、`output`、`columns`（qseqid、sseqid、pident、qlen、slen、length、mismatch、gapopen、qstart、qend、sstart、send、evalue、bitscore）、`hits`、`rows`（每个命中一个字典，最多 100000 个）和 `log_tail`。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import gget_diamond_align
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("proteins.fa").write_text(">p1\nMKTAYIAKQRQISFVKSHFSRQ\n")
    try:
        gget_diamond_align(
            "proteins.fa", "proteins.fa", "hits.tsv",
            diamond_path="no_such_diamond", verbose=False,
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split("] ", 1)[-1].split(" Hint:")[0])
```

- **示例输出：**

```text
BackendUnavailableError diamond_path does not exist: no_such_diamond
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_diamond_align

result = gget_diamond_align("query.fa", "reference.fa", "hits.tsv", sensitivity="sensitive")
print(result["hits"], result["rows"][0]["sseqid"])
```

## 33) `dnakit.databases.gget_genbank_metadata()`

- **作用：** 通过 NCBI E-utilities 获取一组登录号的 GenBank 记录，并写出 XML 与元数据 CSV。
- **处理方法：** 校验登录号，并要求 `xml_path`、`csv_path` 和可选的 `failed_log_path` 是已存在文件夹中的新文件。然后用批大小、请求间隔和可选的 NCBI API key 调用 `gget.gget_virus.fetch_genbank_metadata`；API key 只传给 gget，不会被返回。
- **API：** `gget_genbank_metadata(accessions: 'Sequence[str]', *, xml_path: 'str', csv_path: 'str', batch_size: 'int' = 200, delay: 'float' = 0.5, failed_log_path: 'str | None' = None, api_key: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。NCBI API key 为可选项，可提高请求速率。示例说明：示例传入空的登录号列表，它在检查任何文件路径之前就被拒绝。
- **输入：** `accessions`：1 至 10000 个标识符的列表；`xml_path`、`csv_path`：新的输出文件；`batch_size`：1 至 1000（默认 200）；`delay`：请求间隔秒数，0 至 60（默认 0.5）；`failed_log_path`：可选的新文件；`api_key`：可选的 NCBI key。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`accessions`（数量）、`xml_path`、`csv_path` 和 `summary`（gget 返回的值）。
- **示例代码：**

```python
from dnakit.databases import gget_genbank_metadata
from dnakit.exceptions import ConfigurationError

try:
    gget_genbank_metadata([], xml_path="records.xml", csv_path="records.csv", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError accessions must hold between 1 and 10000 identifiers.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_genbank_metadata

result = gget_genbank_metadata(
    ["NC_045512.2"], xml_path="records.xml", csv_path="records.csv"
)
print(result["summary"])
```

## 34) `dnakit.databases.gget_muscle_align()`

- **作用：** 用 MUSCLE 5 比对核酸或氨基酸序列，并写出比对后的 FASTA 文件。
- **处理方法：** 检查 `fasta` 存在且 `output_path` 是新文件，查找 `muscle` 程序（来自 `muscle_path` 或 PATH；找不到时抛出带安装提示的 `BackendUnavailableError`），运行 `muscle -version` 并要求为 5 版（否则抛出 `BackendVersionError`）。然后运行 `muscle -align <fasta> -output <output_path>`（`super5` 为 True 时用 `-super5`），与 gget 使用的命令相同，并汇总比对后的 FASTA（记录数、比对长度、标识符）。
- **API：** `gget_muscle_align(fasta: 'str', output_path: 'str', *, super5: 'bool' = False, muscle_path: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要外部程序 MUSCLE 5 位于 PATH 中，或在 `muscle_path` 中给出其路径（gget 在 `gget/bins/Linux/muscle` 下附带 Linux 版本）。不会导入 gget 包本身。此环境的 PATH 中没有 MUSCLE，所以示例展示程序缺失错误。示例说明：示例把 `muscle_path` 指向不存在的文件，因此无论是否安装 MUSCLE 都会以 `BackendUnavailableError` 失败。
- **输入：** `fasta`：已存在的 FASTA 文件；`output_path`：新的文件路径，其所在文件夹须存在；`super5`：对大规模输入使用 Super5 算法（默认 False）；`muscle_path`：可选的程序路径。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`program`（名称、版本、路径）、`command`、`output`、`records`、`aligned_length`（序列长度不一致时为 `None`）、`ids` 和 `log_tail`。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import gget_muscle_align
from dnakit.exceptions import BackendUnavailableError

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    Path("seqs.fa").write_text(">a\nACGTACGT\n>b\nACGTACGA\n")
    try:
        gget_muscle_align(
            "seqs.fa", "aligned.fa", muscle_path="no_such_muscle", verbose=False
        )
    except BackendUnavailableError as error:
        print(type(error).__name__, str(error).split("] ", 1)[-1].split(" Hint:")[0])
```

- **示例输出：**

```text
BackendUnavailableError muscle_path does not exist: no_such_muscle
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_muscle_align

result = gget_muscle_align("seqs.fa", "aligned.fa")
print(result["records"], result["aligned_length"], result["ids"])
```

## 35) `dnakit.databases.gget_ref()`

- **作用：** 返回某物种参考文件的 Ensembl FTP 链接，但不下载文件。
- **处理方法：** 校验物种（小写字母、数字和下划线）和 `which`，然后调用 `gget.ref(species, which, release, ftp, save=False)`。gget 在 Ensembl FTP 站点上查找链接并返回，不保存任何文件。`ftp=True` 时 gget 只返回链接。
- **API：** `gget_ref(species: 'str', *, which: 'Sequence[str] | str' = 'all', release: 'int | None' = None, ftp: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入含大写字母和空格的物种名，它在导入 gget 之前就未通过格式检查。
- **输入：** `species`：形如 `homo_sapiens` 的 `genus_species`，或 `human` 这类简称；`which`：`all`，或 `gtf`、`cdna`、`dna`、`cds`、`cdrna`、`pep` 的子集（字符串或列表）；`release`：可选的 Ensembl 发布号，1 至 1000（默认最新）；`ftp`：只返回链接（默认 False）。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`species`、`which`（列表）、`release`、`ftp_only` 和 `links`（gget 返回的值）。
- **示例代码：**

```python
from dnakit.databases import gget_ref
from dnakit.exceptions import ConfigurationError

try:
    gget_ref("Homo Sapiens", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError species has an unsupported format: 'Homo Sapiens'.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_ref

result = gget_ref("homo_sapiens", which=["gtf", "dna"], ftp=True)
print(result["links"])
```

## 36) `dnakit.databases.gget_sars_cov2_download()`

- **作用：** 把 SARS-CoV-2 序列从 NCBI 缓存数据包下载到一个文件夹。
- **处理方法：** 必要时创建 `outdir`，然后用可选筛选条件（宿主、仅完整基因组、是否注释、谱系、登录号以及是否使用登录号）调用 `gget.gget_virus.download_sars_cov2_optimized`。gget 返回 ZIP 路径、NCBI 能够应用的筛选条件，以及无法应用的条件。
- **API：** `gget_sars_cov2_download(*, outdir: 'str', host: 'str | None' = None, complete_only: 'bool | None' = None, annotated: 'bool | None' = None, lineage: 'str | None' = None, accession: 'str | None' = None, use_accession: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入空的文件夹路径，它在导入 gget 或发出任何请求之前就被拒绝。
- **输入：** `outdir`：输出文件夹路径（缺失时创建）；`host`、`lineage`、`accession`：可选字符串；`complete_only`、`annotated`：可选布尔值；`use_accession`：布尔值（默认 False）。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`outdir`、`zip_path`、`applied_filters` 和 `missing_filters`。
- **示例代码：**

```python
from dnakit.databases import gget_sars_cov2_download
from dnakit.exceptions import ConfigurationError

try:
    gget_sars_cov2_download(outdir="", lineage="B.1.1.7", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError outdir must be a folder path.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_sars_cov2_download

result = gget_sars_cov2_download(outdir="sars_cov2", lineage="B.1.1.7", complete_only=True)
print(result["zip_path"], result["missing_filters"])
```

## 37) `dnakit.databases.gget_virus()`

- **作用：** 把病毒基因组数据集（序列和元数据）从 NCBI Virus 下载到一个文件夹。
- **处理方法：** 校验分类单元名称（`is_accession` 为 True 时除外）、序列长度上下限和日期，创建 `outfolder`，并原样带着筛选条件调用 `gget.virus(...)`；未设置的筛选条件保持为 `None`。gget 返回的值被转换成普通 Python 数据，作为 `result` 返回。
- **API：** `gget_virus(virus: 'str', *, outfolder: 'str', is_accession: 'bool' = False, host: 'str | None' = None, min_seq_length: 'int | None' = None, max_seq_length: 'int | None' = None, nuc_completeness: 'str | None' = None, annotated: 'bool | None' = None, geographic_location: 'str | None' = None, min_release_date: 'str | None' = None, max_release_date: 'str | None' = None, lineage: 'str | None' = None, genbank_metadata: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入含不允许字符的分类单元名称，它在创建任何文件夹之前就失败。
- **输入：** `virus`：分类单元名称，例如 `Zika virus`（字母、数字和常见标点，最长 128 个字符），或在 `is_accession=True` 时为登录号列表；`outfolder`：输出文件夹；筛选条件 `host`、`nuc_completeness`、`annotated`、`geographic_location`、`lineage`；`min_seq_length` 与 `max_seq_length`（非负整数，下限不得大于上限）；`min_release_date` 与 `max_release_date`（`YYYY-MM-DD`）；`genbank_metadata`：同时获取 GenBank 元数据。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`virus`、`outfolder` 和 `result`（gget 返回的值）。
- **示例代码：**

```python
from dnakit.databases import gget_virus
from dnakit.exceptions import ConfigurationError

try:
    gget_virus("<bad>", outfolder="virus_data", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError virus must be a taxon name (letters, digits and common punctuation).
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_virus

result = gget_virus(
    "Zika virus",
    outfolder="virus_data",
    host="human",
    min_seq_length=10000,
    min_release_date="2020-01-01",
)
print(result["outfolder"])
```

## 38) `dnakit.databases.gget_virus_metadata()`

- **作用：** 从 NCBI Datasets API 获取病毒基因组元数据，但不下载序列。
- **处理方法：** 校验名称和可选日期，然后用宿主、地点、注释、完整性和 RefSeq 筛选条件调用 `gget.gget_virus.fetch_virus_metadata`。返回的记录被转换成普通数据；gget 返回列表时 `count` 为其数量。
- **API：** `gget_virus_metadata(virus: 'str', *, accession: 'bool' = False, host: 'str | None' = None, geographic_location: 'str | None' = None, annotated: 'bool | None' = None, complete_only: 'bool | None' = None, min_release_date: 'str | None' = None, refseq_only: 'bool | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入空的病毒名称，它在导入 gget 或发出任何请求之前就被拒绝。
- **输入：** `virus`：非空名称（最长 256 个字符）；`accession`：把 `virus` 当作登录号（布尔值，默认 False）；`host`、`geographic_location`：可选字符串；`annotated`、`complete_only`、`refseq_only`：可选布尔值；`min_release_date`：可选的 `YYYY-MM-DD`。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`virus`、`records` 和 `count`。
- **示例代码：**

```python
from dnakit.databases import gget_virus_metadata
from dnakit.exceptions import ConfigurationError

try:
    gget_virus_metadata("", verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError virus must be a non-empty string of at most 256 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_virus_metadata

result = gget_virus_metadata("Zika virus", host="human", complete_only=True)
print(result["count"])
```

## 39) `dnakit.databases.gget_virus_metadata_chunked()`

- **作用：** 按发布日期分块获取超大病毒数据集的元数据。
- **处理方法：** 校验名称和日期，然后调用 `gget.gget_virus.fetch_virus_metadata_chunked`，它把查询拆分为若干发布日期区间以保持每次请求较小。返回的记录被转换成普通数据；gget 返回列表时 `count` 为其数量。
- **API：** `gget_virus_metadata_chunked(virus: 'str', *, accession: 'bool' = False, host: 'str | None' = None, geographic_location: 'str | None' = None, annotated: 'bool | None' = None, complete_only: 'bool' = False, min_release_date: 'str | None' = None, max_release_date: 'str | None' = None, refseq_only: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。示例说明：示例传入数字而不是名称，它在导入 gget 或发出任何请求之前就被拒绝。
- **输入：** `virus`：非空名称（最长 256 个字符）；`accession`：布尔值（默认 False）；`host`、`geographic_location`：可选字符串；`annotated`：可选布尔值；`complete_only`、`refseq_only`：布尔值（默认 False；不接受 `None`）；`min_release_date`、`max_release_date`：可选的 `YYYY-MM-DD`。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`virus`、`records` 和 `count`。
- **示例代码：**

```python
from dnakit.databases import gget_virus_metadata_chunked
from dnakit.exceptions import ConfigurationError

try:
    gget_virus_metadata_chunked(123, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError virus must be a non-empty string of at most 256 characters.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_virus_metadata_chunked

result = gget_virus_metadata_chunked(
    "Severe acute respiratory syndrome coronavirus 2",
    min_release_date="2023-01-01",
    max_release_date="2023-03-31",
)
print(result["count"])
```

## 40) `dnakit.databases.gget_virus_sequences_by_accession()`

- **作用：** 按批次从 NCBI 下载一组登录号对应的病毒基因组序列。
- **处理方法：** 校验登录号，创建 `outdir`，并用批大小和可选的 NCBI API key（只传给 gget，不会被返回）调用 `gget.gget_virus.download_sequences_by_accessions`。gget 返回其写出的 FASTA 文件路径。
- **API：** `gget_virus_sequences_by_accession(accessions: 'Sequence[str]', *, outdir: 'str', batch_size: 'int' = 200, api_key: 'str | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和互联网访问。NCBI API key 为可选项，可提高请求速率。示例说明：示例包含带空格的登录号，它在创建任何文件夹之前就未通过标识符检查。
- **输入：** `accessions`：1 至 10000 个标识符的列表（字母、数字和 `._:/+-`，每个最长 64 个字符）；`outdir`：输出文件夹；`batch_size`：1 至 1000（默认 200）；`api_key`：可选的 NCBI key。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`accessions`（数量）、`outdir` 和 `fasta_path`。
- **示例代码：**

```python
from dnakit.databases import gget_virus_sequences_by_accession
from dnakit.exceptions import ConfigurationError

try:
    gget_virus_sequences_by_accession(
        ["NC_045512.2", "bad accession"], outdir="virus_fasta", verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError accessions[] has an unsupported format: 'bad accession'.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_virus_sequences_by_accession

result = gget_virus_sequences_by_accession(["NC_045512.2"], outdir="virus_fasta")
print(result["fasta_path"])
```

## 41) `dnakit.databases.gget_cellxgene_query()`

- **作用：** 按基因和细胞元数据筛选条件查询 CELLxGENE Census 表达数据，并把结果写到新文件。
- **处理方法：** 先检查参数：`output_path` 必须是新文件，且以 `.h5ad` 结尾（`meta_only` 为 True 时以 `.csv` 结尾），物种必须是五种之一，基因、列名和五类标签筛选条件必须符合允许的格式。然后检查能否导入 `cellxgene_census`，并以 `out=output_path` 调用 `gget.gget_cellxgene.cellxgene(...)`。gget 写出表达量 AnnData（`.h5ad`），`meta_only` 时写出细胞元数据表（`.csv`）；gget 没有返回数据时抛出 `BackendExecutionError`。
- **API：** `gget_cellxgene_query(output_path: 'str', *, species: 'str' = 'homo_sapiens', genes: 'list[str] | None' = None, ensembl: 'bool' = False, tissue: 'list[str] | None' = None, cell_type: 'list[str] | None' = None, disease: 'list[str] | None' = None, sex: 'list[str] | None' = None, assay: 'list[str] | None' = None, column_names: 'list[str] | None' = None, is_primary_data: 'bool' = True, meta_only: 'bool' = False, census_version: 'str' = 'stable', verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 gget（extra `dnakit[bio-gget]`）和 cellxgene-census（extra `dnakit[bio-cellxgene]`），以及对 Census 的互联网访问。此环境未安装 cellxgene-census，所以示例展示在导入该包之前就发生的输入检查。示例说明：示例请求表达量查询（`meta_only` 为 False）却给出 `.csv` 输出；扩展名检查在导入 Census 包之前就失败，因此输出与该包是否安装无关。
- **输入：** `output_path`：新的文件路径；`species`：`homo_sapiens`、`mus_musculus`、`macaca_mulatta`、`callithrix_jacchus` 或 `pan_troglodytes`；`genes`：最多 500 个符号（`ensembl` 为 True 时为 Ensembl ID）；`tissue`、`cell_type`、`disease`、`sex`、`assay`：最多 50 个 Census 标签的列表；`column_names`；`is_primary_data`（默认 True）；`meta_only`（默认 False）；`census_version`（默认 `stable`）。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`species`、`census_version`、`genes`、`ensembl`、`filters`（仅包含已设置的筛选条件）、`is_primary_data`、`meta_only`、`output_path`、`counts`（元数据为 `rows`，否则为 `cells` 和 `genes`）、`obs_columns` 和 `var_names`（各最多 200 项）。
- **示例代码：**

```python
from dnakit.databases import gget_cellxgene_query
from dnakit.exceptions import ConfigurationError

try:
    gget_cellxgene_query("expression.csv", genes=["TP53"], tissue=["lung"], verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError output_path must end in .h5ad for this query.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_cellxgene_query

result = gget_cellxgene_query(
    "lung_tp53.h5ad", genes=["TP53"], tissue=["lung"], cell_type=["macrophage"]
)
print(result["counts"])
```

## 42) `dnakit.databases.gget_cosmic_download_reference()`

- **作用：** 用用户自己的账号下载一个 COSMIC 参考数据包并解压。
- **处理方法：** 检查链接是 `cancer.sanger.ac.uk` 上的 https 地址且其 `path` 参数指向 `.tar` 文件，并检查邮箱与密码。第一次请求以 HTTP 基本认证发送账号，得到带签名的下载地址；随后在不带凭据的情况下从该地址取得数据包，保存为 `<output_folder>/<name>.tar`，经成员检查（拒绝绝对路径、`..` 与链接）后解压到 `<output_folder>/<name>/`。`unzip` 为 True 时，每个解压出的 `.gz` 文件还会就地解压缩。密码不会被保存或返回。传输只使用 Python 标准库，不导入 gget。
- **API：** `gget_cosmic_download_reference(download_link: 'str', output_folder: 'str', *, email: 'str', password: 'str', unzip: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 普通 pip 安装后即可使用（不需要 gget），但需要互联网访问，以及有权下载该文件的 COSMIC 账号（邮箱和密码）。示例说明：示例使用另一个主机上的链接；它在发出任何请求之前就被拒绝，所以示例不会发送任何凭据。
- **输入：** `download_link`：`cancer.sanger.ac.uk` 上的 https 链接，其 `path=` 参数以 `.tar` 结尾；`output_folder`：目标文件夹（缺失时创建）；`email`、`password`：账号（仅限关键字参数）；`unzip`：解压 `.gz` 文件（默认 False）。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`download_link`、`tar_path`、`extracted_folder`、`files`（最多 200 个解压出的路径）和 `unzipped`（已解压缩文件的路径）。
- **示例代码：**

```python
from dnakit.databases import gget_cosmic_download_reference
from dnakit.exceptions import ConfigurationError

try:
    gget_cosmic_download_reference(
        "https://example.org/data.tar",
        "cosmic",
        email="user@example.org",
        password="not-a-real-password",
        verbose=False,
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split("] ", 1)[-1])
```

- **示例输出：**

```text
ConfigurationError download_link must be an https address on cancer.sanger.ac.uk.
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_cosmic_download_reference

link = (
    "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted"
    "?path=grch38/cosmic/v99/Cosmic_MutantCensus_Tsv_v99_GRCh38.tar&bucket=downloads"
)
result = gget_cosmic_download_reference(
    link, "cosmic", email="you@example.org", password="your-cosmic-password", unzip=True
)
print(result["extracted_folder"])
```

## 43) `dnakit.databases.gget_cosmic_select_reference()`

- **作用：** 为某个基因组版本下载一个 COSMIC 项目发布，或在文件已存在时直接复用。
- **处理方法：** 检查项目名（`cancer`、`cell_line`、`census`、`resistance`、`genome_screen` 或 `targeted_screen`）、GRCh 版本（37 或 38；`cancer` 只有 37）和发布号（1 至 200），构造数据包名与预期文件名，并在 `output_folder` 中查找。若 `<output_folder>/<数据包名>/<文件>` 已存在（`.tsv`，`unzip` 为 False 时为 `.tsv.gz`）且 `overwrite` 为 False，则立即返回，`downloaded` 为 False。否则用该账号调用 `gget_cosmic_download_reference` 并返回文件路径。
- **API：** `gget_cosmic_select_reference(output_folder: 'str', cosmic_project: 'str', grch_version: 'int', cosmic_version: 'int', *, email: 'str', password: 'str', unzip: 'bool' = True, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **运行要求：** 普通 pip 安装后即可使用（不需要 gget）。复用已有文件不需要网络；下载则需要互联网访问和 COSMIC 账号（邮箱和密码）。示例说明：示例创建了函数为 `census` 项目、GRCh38、第 99 版所预期的文件，因此函数直接复用它，不发出任何请求。账号参数只是占位符，也不会被返回。
- **输入：** `output_folder`：文件所在文件夹；`cosmic_project`：项目名；`grch_version`：37 或 38；`cosmic_version`：发布号；`email`、`password`：账号（仅限关键字参数；复用已有文件时不会用到）；`unzip`：默认 True；`overwrite`：文件已存在时仍重新下载（默认 False）。
- **输出：** 字典，含gget 通用结果字段 `operation`、`status`（`completed`）、`backend`（`gget`）、`package`（`gget <版本>`）、`service`（所用服务的名称）和 `accessed`（当天日期）、`cosmic_project`、`grch_version`、`cosmic_version`、`file_path` 和 `downloaded`。复用时还含 `overwritten`（False）；下载后则改为含 `download_link`、`overwritten`、`tar_path` 和 `file_count`。
- **示例代码：**

```python
import os
import tempfile
from pathlib import Path

from dnakit.databases import gget_cosmic_select_reference

with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    release = Path("cosmic") / "Cosmic_MutantCensus_Tsv_v99_GRCh38"
    release.mkdir(parents=True)
    (release / "Cosmic_MutantCensus_v99_GRCh38.tsv").write_text(
        "GENE_SYMBOL\tMUTATION_AA\nTP53\tp.R175H\n"
    )
    result = gget_cosmic_select_reference(
        "cosmic", "census", 38, 99,
        email="user@example.org", password="not-a-real-password", verbose=False,
    )
    print(result["operation"], result["cosmic_project"], result["grch_version"], result["cosmic_version"])
    print(Path(result["file_path"]).name, result["downloaded"], result["overwritten"])
    print("password" in str(result))
```

- **示例输出：**

```text
select_reference census 38 99
Cosmic_MutantCensus_v99_GRCh38.tsv False False
False
```

- **安装依赖后的用法：**

```python
from dnakit.databases import gget_cosmic_select_reference

result = gget_cosmic_select_reference(
    "cosmic", "census", 38, 99, email="you@example.org", password="your-cosmic-password"
)
print(result["file_path"], result["downloaded"])
```
