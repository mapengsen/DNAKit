# fastp 读段处理与 Kraken 2 物种分类

两组封装用户自行安装的程序的命令行函数：[fastp](https://github.com/OpenGene/fastp) 的读段剪切与过滤（包含 `run_fastp` 未暴露的选项，位于 `dnakit.sequencing_qc`），以及 [Kraken 2](https://github.com/DerrickWood/kraken2) 的分类、数据库检查、合并和数据库构建（位于 `dnakit.microbiome`）。安装：`conda install -c bioconda fastp kraken2`。

每个函数用类型化参数构造固定的参数列表来运行一个程序：不经过 shell，没有通用执行器，输出都是新文件，只传入你设置的选项，其余选项保持程序自身默认值。找不到程序时抛出带安装提示的 `BackendUnavailableError`；`fastp_path` 指定 fastp 可执行文件，`kraken2_path` 指定包含 `kraken2`、`kraken2-build`、`k2` 的目录。返回字典包含 `status`（退出码为 0 且所有预期输出存在时才为 `completed`）、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 以及少量解析结果。测试版本：fastp 1.4.0、Kraken 2 2.17.2。

`run_fastp` 已覆盖 JSON/HTML 报告、接头序列、双端接头检测、poly-X 剪切、滑动窗口质量剪切和低质量碱基比例限制；`fastp_process_reads` 增加接头 FASTA、低复杂度与 index 过滤、重叠区纠错、UMI 处理和过表达序列分析（同时也接受常用选项）。`classify_kraken2` 使用默认选项分类；`kraken2_classify_reads` 增加报告、快速模式、多数据库及其他选项。下载与构建函数需要访问 NCBI 或 GTDB 服务器，可能耗时数小时，因此这里的示例用本地序列构建一个小型数据库。

## 示例数据

准备示例为 fastp 写入 300 对双端 reads（其中 100 对插入片段很短，reads 读进 Illumina 接头 `AGATCGGAAGAGC`）和接头 FASTA；为 Kraken 2 写入一个玩具分类树（root、Bacteria、Genusa 属及 alpha、beta 两个种）、两个带分类号标记的 4 kb 随机基因组，以及 45 条 reads（每个基因组 20 条，另加 5 条随机序列）。

- **示例代码：**

```python
import random

rng = random.Random(9)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


rc = str.maketrans("ACGT", "TGCA")
adapter = "AGATCGGAAGAGC"
first, second = [], []
for index in range(300):
    fragment = dna(120)
    read1, read2 = fragment[:100], fragment[::-1].translate(rc)[:100]
    if index < 100:
        insert = fragment[:70]
        read1 = (insert + adapter + "ACGTACGTAC" * 3)[:100]
        read2 = (insert[::-1].translate(rc) + adapter + "ACGTACGTAC" * 3)[:100]
    first.append(f"@p{index}/1\n{read1}\n+\n{'I' * 100}\n")
    second.append(f"@p{index}/2\n{read2}\n+\n{'I' * 100}\n")
open("r1.fq", "w").write("".join(first))
open("r2.fq", "w").write("".join(second))
open("adapters.fa", "w").write(f">illumina\n{adapter}\n")
import os

os.makedirs("db/taxonomy")
open("db/taxonomy/nodes.dmp", "w").write("1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n10\t|\t2\t|\tgenus\t|\n11\t|\t10\t|\tspecies\t|\n12\t|\t10\t|\tspecies\t|\n")
names = ((1, "root"), (2, "Bacteria"), (10, "Genusa"), (11, "Genusa alpha"), (12, "Genusa beta"))
open("db/taxonomy/names.dmp", "w").write("".join(f"{i}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for i, n in names))
genomes = {11: dna(4000), 12: dna(4000)}
for taxid, sequence in genomes.items():
    open(f"genome_{taxid}.fa", "w").write(f">seq{taxid}|kraken:taxid|{taxid}\n{sequence}\n")
reads = "".join(f">a{i}\n{genomes[11][s:s + 150]}\n" for i, s in enumerate(range(0, 3000, 150)))
reads += "".join(f">b{i}\n{genomes[12][s:s + 150]}\n" for i, s in enumerate(range(0, 3000, 150)))
reads += "".join(f">u{i}\n{dna(150)}\n" for i in range(5))
open("reads.fa", "w").write(reads)
print("ready")
```

- **示例输出：**

```text
ready
```

## 1) `dnakit.sequencing_qc.fastp_process_reads()`

- **作用：** 用 fastp 剪切并过滤 FASTQ reads，涵盖接头 FASTA、低复杂度与 index 过滤、重叠区纠错、UMI 处理和过表达序列分析。
- **处理方法：** 运行 `fastp -i ... -o ...`（双端再加 `-I`/`-O`），选项由参数构成：`--adapter_fasta`、`--adapter_sequence`、`--detect_adapter_for_pe`、`--low_complexity_filter` 与 `--complexity_threshold`、`--filter_by_index1/2`、`--correction`、`--umi` 与 `--umi_loc`、`--umi_len`、`--umi_prefix`、`--overrepresentation_analysis`、poly-X 与滑动窗口剪切以及质量过滤。JSON 报告会被解析为 `summary`。
- **API：** `dnakit.sequencing_qc.fastp_process_reads(read1, output_read1, json_report, *, read2=None, output_read2=None, html_report=None, adapter_sequence=None, adapter_sequence_read2=None, adapter_fasta=None, detect_adapter_for_pe=False, disable_adapter_trimming=False, low_complexity_filter=False, complexity_threshold=None, filter_by_index1=None, filter_by_index2=None, index_filter_threshold=None, correction=False, overlap_length_required=None, overlap_diff_limit=None, overlap_diff_percent_limit=None, umi=False, umi_location=None, umi_length=None, umi_prefix=None, umi_skip=None, umi_delimiter=None, overrepresentation_analysis=False, overrepresentation_sampling=None, trim_poly_x=False, cut_front=False, cut_tail=False, cut_right=False, cut_window_size=None, cut_mean_quality=None, qualified_quality_phred=None, unqualified_percent_limit=None, n_base_limit=None, average_quality=None, length_required=None, threads=None, fastp_path=None, timeout_seconds=86400.0)`.
- **输入：** read 1（双端时加 read 2 与 `output_read2`）、新的输出路径和 `json_report`，以及上述可选参数和 `html_report`、`threads`、`timeout_seconds`。
- **输出：** 通用结果键，加上 `summary`（过滤前后统计、过滤结果计数、重复率、接头切除、过表达序列）。
- **示例代码：**

```python
from dnakit.sequencing_qc import fastp_p14_cli_functions as fp

result = fp.fastp_process_reads("r1.fq", "clean_1.fq", "report.json", read2="r2.fq", output_read2="clean_2.fq", html_report="report.html",
                                adapter_fasta="adapters.fa", threads=2, verbose=False)
summary = result["summary"]
lengths = {len(line) for line in open("clean_1.fq").read().split("\n")[1::4] if line}
print(result["status"], summary["before_filtering"]["total_reads"], summary["after_filtering"]["total_reads"], summary["adapter_cutting"]["adapter_trimmed_reads"])
print(min(lengths), max(lengths))
```

- **示例输出：**

```text
completed 600 600 201
70 100
```

## 2) `dnakit.microbiome.kraken2_classify_reads()`

- **作用：** 把 reads 对一个或多个 Kraken 2 数据库分类，并写出各分类单元的 clade 计数报告。
- **处理方法：** 运行 `kraken2 --db <库> --output <文件>`（按需加 `--report`、`--quick`、`--confidence`、`--minimum-hit-groups`、`--paired`、`--report-minimizer-data` 等）。多个数据库、给出 `merge_policy` 或 `interleaved` 时改用 `k2 classify`，默认按最低共同祖先合并各库的结果。蛋白数据库用同一调用对 DNA reads 做六框翻译搜索。已处理/已分类计数从日志读取，报告被解析为 `report_rows`。
- **API：** `dnakit.microbiome.kraken2_classify_reads(reads, databases, output, *, report=None, paired=False, interleaved=False, quick=False, confidence=None, minimum_hit_groups=None, minimum_base_quality=None, merge_policy=None, threads=None, use_names=False, memory_mapping=False, report_minimizer_data=False, report_zero_counts=False, use_mpa_style=False, classified_out=None, unclassified_out=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **输入：** reads 文件列表、数据库目录（或其列表）、新的输出路径，以及可选的 `report`、`paired`、`interleaved`、`quick`、`confidence`、`minimum_hit_groups`、`minimum_base_quality`、`merge_policy`、`threads`、`use_names`、`memory_mapping`、`report_minimizer_data`、`report_zero_counts`、`use_mpa_style`、`classified_out`、`unclassified_out`。
- **输出：** 通用结果键，加上 `databases`、`sequences_processed`、`sequences_classified`、`sequences_unclassified`、`percent_classified` 和 `report_rows`（percent、clade_reads、taxon_reads、rank、taxid、name）。
- **示例代码：**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

added = k2.kraken2_add_to_library("db", ["genome_11.fa", "genome_12.fa"], no_masking=True, verbose=False)
built = k2.kraken2_build_database("db", kmer_length=31, minimizer_length=15, minimizer_spaces=3, threads=2, verbose=False)
result = k2.kraken2_classify_reads(["reads.fa"], "db", "classified.txt", report="report.txt", verbose=False)
print(added["status"], built["status"], result["status"])
print(result["sequences_processed"], result["sequences_classified"], result["sequences_unclassified"])
for row in result["report_rows"]:
    print(row["rank"], row["taxid"], row["clade_reads"], row["taxon_reads"], row["name"])
```

- **示例输出：**

```text
completed completed completed
45 40 5
U 0 5 5 unclassified
R 1 40 0 root
D 2 40 0 Bacteria
G 10 40 0 Genusa
S 11 20 20 Genusa alpha
S 12 20 20 Genusa beta
```

## 3) `dnakit.microbiome.kraken2_add_to_library()`

- **作用：** 把自己的 FASTA 序列加入 Kraken 2 数据库的 library。
- **处理方法：** `backend='kraken2-build'` 时对每个文件运行一次 `kraken2-build --add-to-library <文件> --db <库>`，`backend='k2'` 时运行一次 `k2 add-to-library --file ...`。每条序列的标题中需要分类号（`>名称|kraken:taxid|<id>`），或其登录号列在 `seqid2taxid.map` 中。除非 `no_masking`，否则会屏蔽低复杂度序列。
- **API：** `dnakit.microbiome.kraken2_add_to_library(database, fasta_files, *, backend='kraken2-build', protein=False, no_masking=False, threads=None, masker_threads=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **输入：** 数据库目录（不存在则创建）、FASTA 文件列表，以及可选的 `backend`、`protein`、`no_masking`、`threads`、`masker_threads`（仅 k2）。
- **输出：** 通用结果键、`runs`（`kraken2-build` 时每个文件一个结果）和 `outputs`（`library` 目录）。
- **示例代码：**

```python
import os

from dnakit.microbiome import kraken2_p14_cli_functions as k2

os.makedirs("db_k2/taxonomy")
for name in ("nodes.dmp", "names.dmp"):
    open(f"db_k2/taxonomy/{name}", "w").write(open(f"db/taxonomy/{name}").read())
result = k2.kraken2_add_to_library("db_k2", ["genome_11.fa"], backend="k2", no_masking=True, verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]], sorted(os.listdir("db_k2/library")))
```

- **示例输出：**

```text
completed ['library'] ['added']
```

## 4) `dnakit.microbiome.kraken2_build_database()`

- **作用：** 根据 library 构建 Kraken 2 数据库的 minimizer 哈希表和分类索引，也可构建用于翻译搜索的蛋白数据库。
- **处理方法：** 运行 `kraken2-build --build --db <库>`（或 `k2 build --db <库>`），按需加 k-mer、minimizer、spaces、负载因子、大小、线程、`--fast-build` 和 `--protein` 选项。目录中需有 `taxonomy/` 和 `library/`。`protein=True` 时，该库通过六框翻译搜索对 DNA reads 分类。
- **API：** `dnakit.microbiome.kraken2_build_database(database, *, backend='kraken2-build', kmer_length=None, minimizer_length=None, minimizer_spaces=None, load_factor=None, max_db_size=None, threads=None, fast_build=False, protein=False, parts=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **输入：** 数据库目录，以及可选的 `backend`、`kmer_length`（默认 35，蛋白 15）、`minimizer_length`（31，12）、`minimizer_spaces`（7，0）、`load_factor`（0.7）、`max_db_size`、`threads`、`fast_build`、`protein`、`parts`（仅 k2）。
- **输出：** 通用结果键；`outputs` 列出 `hash.k2d`、`opts.k2d`、`taxo.k2d`。
- **示例代码：**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_build_database("db_k2", backend="k2", kmer_length=31, minimizer_length=15, minimizer_spaces=3, threads=2, verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **示例输出：**

```text
completed ['hash.k2d', 'opts.k2d', 'taxo.k2d']
```

## 5) `dnakit.microbiome.kraken2_inspect_database()`

- **作用：** 汇总 Kraken 2 数据库及各分类单元的 minimizer 计数。
- **处理方法：** 运行 `k2 inspect --db <库> --output <文件>`，其报告与 `kraken2-inspect` 相同：头部行（k-mer 与 minimizer 长度、间隔掩码、分类节点数、表大小），除非 `skip_counts`，其后每个分类单元一行。
- **API：** `dnakit.microbiome.kraken2_inspect_database(database, output, *, skip_counts=False, use_mpa_style=False, report_zero_counts=False, memory_mapping=False, kraken2_path=None, timeout_seconds=3600.0)`.
- **输入：** 已构建的数据库目录、新的输出路径，以及可选的 `skip_counts`、`use_mpa_style`、`report_zero_counts`、`memory_mapping`。
- **输出：** 通用结果键，加上 `header`（`#` 行）和 `rows`（percent、clade_minimizers、taxon_minimizers、rank、taxid、name）。
- **示例代码：**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_inspect_database("db", "inspect.txt", verbose=False)
print(result["status"], result["header"][1])
for row in result["rows"]:
    print(row["rank"], row["taxid"], row["name"])
```

- **示例输出：**

```text
completed Database options: nucleotide db, k = 31, l = 15
R 1 root
D 2 Bacteria
G 10 Genusa
S 12 Genusa beta
S 11 Genusa alpha
```

## 6) `dnakit.microbiome.kraken2_merge_classifications()`

- **作用：** 合并同一批 reads 对不同数据库的逐条分类结果，冲突时取最低共同祖先。
- **处理方法：** 运行 `k2 merge-classifications --db <库1>,<库2> --output <文件> <文件...>`（按需加 `--report`）。
- **API：** `dnakit.microbiome.kraken2_merge_classifications(classification_files, databases, output, *, report=None, use_names=False, use_mpa_style=False, report_zero_counts=False, threads=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **输入：** 至少两个分类结果文件、产生它们的数据库、新的输出路径，以及可选的 `report`、`use_names`、`use_mpa_style`、`report_zero_counts`、`threads`。
- **输出：** 通用结果键；写出合并后的分类结果（及报告）。
- **示例代码：**

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

k2.kraken2_classify_reads(["reads.fa"], "db", "quick.txt", quick=True, verbose=False)
result = k2.kraken2_merge_classifications(["classified.txt", "quick.txt"], ["db", "db"], "merged.txt", report="merged_report.txt", verbose=False)
print(result["status"], len(open("merged.txt").read().splitlines()))
```

- **示例输出：**

```text
completed 45
```

## 7) `dnakit.microbiome.kraken2_download_library()`

- **作用：** 把 NCBI 基因组序列下载到 Kraken 2 数据库的 library 目录。
- **处理方法：** 运行 `kraken2-build --download-library <类型> --db <库>` 或 `k2 download-library --library <名称|taxid|project|accession>`，按需加屏蔽、FTP、组装来源、组装水平、注释、续传和 BLAST 分卷选项。构建前需另行下载分类树。
- **API：** `dnakit.microbiome.kraken2_download_library(database, library, *, backend='kraken2-build', protein=False, no_masking=False, use_ftp=False, threads=None, assembly_source=None, assembly_levels=None, has_annotation=False, resume=False, blast_volumes=None, masker_threads=None, kraken2_path=None, timeout_seconds=86400.0)`.
- **输入：** 数据库目录（不存在则创建）、`library`（archaea、bacteria、plasmid、viral、human、fungi、plant、protozoa、nr、nt、UniVec、UniVec_Core；`backend='k2'` 时还可用 taxid、project 或 accession），以及可选的 `backend`、`protein`、`no_masking`、`use_ftp`、`threads`，和仅 k2 支持的 `assembly_source`、`assembly_levels`、`has_annotation`、`resume`、`blast_volumes`、`masker_threads`。
- **输出：** 通用结果键；`outputs` 为 `library` 目录。
- **示例调用（此处不运行）：** 该调用需联网从 NCBI 下载，因此不在文档测试中运行；已用 `UniVec_Core` 手工验证（`kraken2-build --use-ftp` 和 `k2 download-library` 两种方式），设置 `DNAKIT_LIVE_TESTS=1` 时单元测试会运行它。

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_download_library("univec_db", "UniVec_Core", use_ftp=True, no_masking=True)
print(result["status"], result["outputs"])
```

## 8) `dnakit.microbiome.kraken2_build_standard_database()`

- **作用：** 下载并构建 Kraken 2 标准数据库。
- **处理方法：** 运行 `kraken2-build --standard --db <库>` 或 `k2 build --standard --db <库>`。标准库包含 RefSeq 的古菌、细菌、病毒、质粒、人类和 UniVec_Core 序列及 NCBI 分类树；下载量数十 GB，构建需要数十 GB 内存和数小时。
- **API：** `dnakit.microbiome.kraken2_build_standard_database(database, *, backend='kraken2-build', no_masking=False, use_ftp=False, kmer_length=None, minimizer_length=None, minimizer_spaces=None, load_factor=None, max_db_size=None, threads=None, fast_build=False, protein=False, kraken2_path=None, timeout_seconds=86400.0)`.
- **输入：** 数据库目录，以及可选的 `backend`、`no_masking`、`use_ftp`、`kmer_length`、`minimizer_length`、`minimizer_spaces`、`load_factor`、`max_db_size`、`threads`、`fast_build`、`protein`、`timeout_seconds`。
- **输出：** 通用结果键；`outputs` 列出 `hash.k2d`、`opts.k2d`、`taxo.k2d`。
- **示例调用（此处不运行）：** 该调用下载数十 GB 并构建数小时，因此不在文档测试中运行。

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

result = k2.kraken2_build_standard_database("standard_db", threads=16, use_ftp=True, timeout_seconds=86400)
print(result["status"])
```

## 9) `dnakit.microbiome.kraken2_build_special_database()`

- **作用：** 下载并构建特殊参考数据库：16S rRNA 数据库 Greengenes、SILVA、RDP，或基因组分类数据库（GTDB）。
- **处理方法：** 运行 `kraken2-build --special <类型> --db <库>` 或 `k2 build --special <类型> --db <库>`。`special='gtdb'` 及 `gtdb_*` 选项（`--gtdb-files`、`--gtdb-use-ncbi-taxonomy`、`--gtdb-representatives-only`、`--gtdb-split-database`、`--gtdb-server`）需要 `backend='k2'`。RDP 的下载地址上游已失效。
- **API：** `dnakit.microbiome.kraken2_build_special_database(database, special, *, backend='kraken2-build', kmer_length=None, minimizer_length=None, minimizer_spaces=None, load_factor=None, max_db_size=None, threads=None, fast_build=False, gtdb_files=None, gtdb_use_ncbi_taxonomy=False, gtdb_representatives_only=False, gtdb_split_database=False, gtdb_server=None, no_masking=False, kraken2_path=None, timeout_seconds=86400.0)`.
- **输入：** 数据库目录、`special`（`greengenes`、`silva`、`rdp`、`gtdb`），以及可选的 `backend`、`kmer_length`、`minimizer_length`、`minimizer_spaces`、`load_factor`、`max_db_size`、`threads`、`fast_build`、`gtdb_files`、`gtdb_use_ncbi_taxonomy`、`gtdb_representatives_only`、`gtdb_split_database`、`gtdb_server`、`no_masking`、`timeout_seconds`。
- **输出：** 通用结果键；`outputs` 列出 `hash.k2d`、`opts.k2d`、`taxo.k2d`。
- **示例调用（此处不运行）：** 这些调用需联网下载 16S 或 GTDB 数据并构建很久，因此不在文档测试中运行。

```python
from dnakit.microbiome import kraken2_p14_cli_functions as k2

silva = k2.kraken2_build_special_database("silva_db", "silva", threads=8)
gtdb = k2.kraken2_build_special_database("gtdb_db", "gtdb", backend="k2", gtdb_representatives_only=True, threads=16)
print(silva["status"], gtdb["status"])
```
