# 比对与搜索程序：BLAST+、BamTools、HISAT2

`dnakit.alignment` 中的这组函数调用用户自行安装的 [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/)、[BamTools](https://github.com/pezmaster31/bamtools) 和 [HISAT2](https://daehwankimlab.github.io/hisat2/)（例如 `conda install -c bioconda blast bamtools hisat2`）。每个公开函数用类型化参数构造固定的参数列表来运行一个程序或子命令：不经过 shell，没有通用执行器，输出都是新文件，只传入你设置的选项，其余选项保持程序自身默认值。找不到程序时抛出带安装提示的 `BackendUnavailableError`；`blast_path`、`bamtools_path`、`hisat2_path` 用来指定程序所在目录（或程序本身），代替 PATH。

返回字典包含 `status`（退出码为 0 且所有预期输出存在时才为 `completed`）、`backend`、`version`、`command`、`return_code`、`outputs`、`log_tail` 以及少量解析结果。测试版本：BLAST+ 2.17.0、BamTools 2.5.3、HISAT2 2.2.3。把多个 BLAST+ 任务或子命令合并为一个参数的函数会写明覆盖哪些取值。

## 示例数据

下面的准备示例用 6 个随机的 900 bp 基因建立小型 BLAST 核酸库 `genes`，从 `cds2` 截取 200 bp 作为查询序列，从 `cds1` 截取 22 bp 作为引物，生成含二核苷酸重复的 40 kb 基因组，以及含 60 条 reads 的 BAM（40 条比对到两条参考序列，20 条未比对）。后面的示例在同一目录中使用这些文件。

- **示例代码：**

```python
import random
import subprocess

import pysam

rng = random.Random(1)


def dna(length):
    return "".join(rng.choice("ACGT") for _ in range(length))


genes = [dna(900) for _ in range(6)]
with open("genes.fa", "w") as handle:
    for index, sequence in enumerate(genes):
        handle.write(f">cds{index}\n{sequence}\n")
subprocess.run(["makeblastdb", "-in", "genes.fa", "-dbtype", "nucl", "-out", "genes", "-parse_seqids"], check=True, capture_output=True)
open("query.fa", "w").write(f">q0\n{genes[2][100:300]}\n")
open("primer.fa", "w").write(f">primer\n{genes[1][50:72]}\n")
genome = dna(20000) + "AC" * 400 + dna(20000)
open("genome.fa", "w").write(">chr1\n" + "\n".join(genome[i:i + 70] for i in range(0, len(genome), 70)) + "\n")
header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 5000}, {"SN": "chr2", "LN": 3000}]}
with pysam.AlignmentFile("raw.bam", "wb", header=header) as out:
    for index in range(60):
        read = pysam.AlignedSegment(out.header)
        read.query_name = f"r{index}"
        read.query_sequence = dna(50)
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        if index < 40:
            read.reference_id = index % 2
            read.reference_start = 100 + 10 * index
            read.cigar = [(0, 50)]
            read.flag = 16 if index % 4 == 1 else 0
            read.mapping_quality = 60 if index % 3 else 10
            read.set_tag("XX", "a" if index % 2 else "b")
        else:
            read.flag = 4
            read.reference_id = -1
        out.write(read)
pysam.sort("-o", "reads.bam", "raw.bam")
pysam.index("reads.bam")
print("ready")
```

- **示例输出：**

```text
ready
```

## 1) `dnakit.alignment.blastn_search()`

- **作用：** 用任意 `blastn` 任务对 BLAST 数据库或 subject FASTA 做核酸-核酸 BLAST 比对。
- **处理方法：** 运行 `blastn`，`-task` 取 `blastn`（经典，词长 11）、`blastn-short`（50 碱基以下的查询，如引物）、`dc-megablast`（不连续 megablast，跨物种）、`megablast`（BLAST+ 默认）或 `rmblastn`（重复元件打分，需要含 RMBlast 的 BLAST+ 构建）。输出格式 6、7、10 会解析为 `rows`。
- **API：** `dnakit.alignment.blastn_search(query, output, *, database=None, subject=None, task='megablast', outfmt=None, evalue=None, word_size=None, max_target_seqs=None, max_hsps=None, perc_identity=None, strand=None, dust=None, gapopen=None, gapextend=None, penalty=None, reward=None, template_type=None, template_length=None, num_threads=None, blast_path=None)`.
- **输入：** 查询 FASTA、新的输出路径、`database`（makeblastdb 生成的前缀）与 `subject`（FASTA）二选一、`task`、`outfmt`（0-18 或如 `6 qseqid sseqid pident` 的字符串），以及可选的 `evalue`、`word_size`、`max_target_seqs`、`max_hsps`、`perc_identity`、`strand`、`dust`、`gapopen`、`gapextend`、`penalty`、`reward`、`template_type`、`template_length`、`num_threads`。
- **输出：** 通用结果键，加上 `task`，表格格式时还有 `rows` / `row_count`。
- **示例代码：**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

for task in ("blastn", "megablast", "dc-megablast"):
    result = blast.blastn_search("query.fa", f"{task}.tsv", database="genes", task=task, outfmt=6, verbose=False)
    top = result["rows"][0]
    print(task, result["status"], top["sseqid"], top["pident"], top["length"])
short = blast.blastn_search("primer.fa", "short.tsv", subject="genes.fa", task="blastn-short", outfmt="6 qseqid sseqid length", verbose=False)
print(short["rows"][0])
```

- **示例输出：**

```text
blastn completed cds2 100.000 200
megablast completed cds2 100.000 200
dc-megablast completed cds2 100.000 200
{'qseqid': 'primer', 'sseqid': 'cds1', 'length': '22'}
```

## 2) `dnakit.alignment.rpstblastn_search()`

- **作用：** 把 DNA 查询序列做六框翻译后，比对到保守结构域 PSSM 数据库。
- **处理方法：** 运行 `rpstblastn`，传入查询序列、RPS-BLAST 数据库前缀（如 NCBI 的 `Cdd`，由 `makeprofiledb` 生成）、遗传密码表和链选项。
- **API：** `dnakit.alignment.rpstblastn_search(query, database, output, *, outfmt=None, evalue=None, max_target_seqs=None, query_gencode=None, strand=None, num_threads=None, blast_path=None)`.
- **输入：** 查询 FASTA、数据库前缀、新的输出路径，以及可选的 `outfmt`、`evalue`、`max_target_seqs`、`query_gencode`、`strand`、`num_threads`。
- **输出：** 通用结果键；表格格式时还有 `rows` / `row_count`。
- **示例调用（此处不运行）：** 该调用需要 NCBI 保守结构域数据库（`Cdd`，用 `update_blastdb.pl` 下载，数 GB），因此不在文档测试中运行。

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.rpstblastn_search("contig.fa", "Cdd", "domains.tsv", outfmt=6, evalue=0.01, query_gencode=1)
print(result["status"], result["row_count"])
```

## 3) `dnakit.alignment.blastdbcmd_fetch_sequences()`

- **作用：** 按标识符从 BLAST 数据库取出序列或其子区间。
- **处理方法：** 运行 `blastdbcmd`，用 `-entry`（逗号连接的标识符，或 `all`）并可加 `-range`、`-strand`；或用 `-entry_batch`（每行一个标识符、可跟区间的文件）。
- **API：** `dnakit.alignment.blastdbcmd_fetch_sequences(database, output, *, entries=None, entry_batch=None, database_type='nucl', sequence_range=None, strand=None, outfmt=None, target_only=False, line_length=None, long_seqids=False, blast_path=None)`.
- **输入：** 数据库前缀、新的输出路径、`entries` 与 `entry_batch` 二选一，以及可选的 `database_type`、`sequence_range`（`10-50`、`10-`）、`strand`、`outfmt`（默认 `%f`，即 FASTA）、`target_only`、`line_length`、`long_seqids`。
- **输出：** 通用结果键；序列写在输出文件中。
- **示例代码：**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

one = blast.blastdbcmd_fetch_sequences("genes", "one.fa", entries=["cds3"], sequence_range="11-40", verbose=False)
print(one["status"], open("one.fa").read().split("\n")[0])
open("ids.txt", "w").write("cds0\ncds4\n")
batch = blast.blastdbcmd_fetch_sequences("genes", "batch.txt", entry_batch="ids.txt", outfmt="%a %l", verbose=False)
print(open("batch.txt").read().split("\n")[:2])
```

- **示例输出：**

```text
completed >cds3:11-40 
['cds0 900', 'cds4 900']
```

## 4) `dnakit.alignment.makembindex_build_index()`

- **作用：** 为核酸 BLAST 数据库或 FASTA 文件建立 MegaBLAST 索引。
- **处理方法：** 运行 `makembindex`。`input_format='blastdb'`（默认）时按新索引格式把 `<prefix>.NN.idx` 写在数据库旁边；`'fasta'` 时按已弃用的旧格式写到 `output_prefix` 下。索引通过 blastn 的 `-use_index` 使用。
- **API：** `dnakit.alignment.makembindex_build_index(input_path, *, input_format='blastdb', output_prefix=None, legacy=None, idmap=False, nmer=None, ws_hint=None, volsize=None, stride=None, blast_path=None)`.
- **输入：** 数据库前缀（或 FASTA 文件），以及可选的 `input_format`、`output_prefix`（仅 fasta）、`legacy`、`idmap`、`nmer`、`ws_hint`、`volsize`、`stride`。
- **输出：** 通用结果键；`outputs` 列出索引分卷。
- **示例代码：**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.makembindex_build_index("genes", verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **示例输出：**

```text
completed ['genes.00.idx']
```

## 5) `dnakit.alignment.windowmasker_count_units()`

- **作用：** 生成 WindowMasker 做重复序列屏蔽所需的 k-mer 单元计数。
- **处理方法：** 对基因组序列运行 `windowmasker -mk_counts`，按所选格式写出计数文件。
- **API：** `dnakit.alignment.windowmasker_count_units(input_path, output, *, unit_length=None, genome_size=None, input_format='fasta', input_is_list=False, check_duplicates=False, memory_mb=None, counts_format='ascii', counts_size_mb=None, metadata=None, blast_path=None)`.
- **输入：** FASTA（或 `input_is_list` 时的 FASTA 文件清单）、新的输出路径，以及可选的 `unit_length`（1-16）、`genome_size`、`input_format`、`check_duplicates`、`memory_mb`、`counts_format`（`ascii`、`binary`、`oascii`、`obinary`）、`counts_size_mb`、`metadata`。
- **输出：** 通用结果键；计数写在输出文件中。
- **示例代码：**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.windowmasker_count_units("genome.fa", "counts.txt", verbose=False)
print(result["status"], open("counts.txt").read().splitlines()[0])
```

- **示例输出：**

```text
completed ##windowmasker-statistics-format-version:ascii 1.1.0
```

## 6) `dnakit.alignment.windowmasker_convert_counts()`

- **作用：** 在 ascii、二进制和优化格式之间转换 WindowMasker 单元计数文件。
- **处理方法：** 运行 `windowmasker -convert`，目标格式由 `-sformat` 指定。优化格式（`oascii`、`obinary`）加载更快。
- **API：** `dnakit.alignment.windowmasker_convert_counts(input_path, output, *, counts_format='ascii', counts_size_mb=None, metadata=None, blast_path=None)`.
- **输入：** 计数文件、新的输出路径、`counts_format`，以及可选的 `counts_size_mb`、`metadata`。
- **输出：** 通用结果键；转换后的计数写在输出文件中。
- **示例代码：**

```python
import os

from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.windowmasker_convert_counts("counts.txt", "counts.obin", counts_format="obinary", verbose=False)
print(result["status"], os.path.getsize("counts.obin") > 0, open("counts.obin", "rb").read(1)[0])
```

- **示例输出：**

```text
completed True 3
```

## 7) `dnakit.alignment.windowmasker_mask_sequences()`

- **作用：** 用 WindowMasker 屏蔽基因组序列中的重复序列。
- **处理方法：** 运行 `windowmasker -ustat <counts>`，输出区间（默认）、软屏蔽 FASTA 或 maskinfo/seqloc 格式之一。留为 None 的阈值沿用计数文件中保存的值。
- **API：** `dnakit.alignment.windowmasker_mask_sequences(input_path, unit_counts, output, *, input_format='fasta', output_format='interval', window=None, t_thres=None, t_extend=None, t_low=None, t_high=None, set_t_low=None, set_t_high=None, dust=False, dust_level=None, parse_seqids=False, blast_path=None)`.
- **输入：** FASTA、`windowmasker_count_units` 生成的计数文件、新的输出路径、`input_format`、`output_format`，以及可选的 `window`、`t_thres`、`t_extend`、`t_low`、`t_high`、`set_t_low`、`set_t_high`、`dust`、`dust_level`、`parse_seqids`。
- **输出：** 通用结果键；屏蔽结果写在输出文件中。
- **示例代码：**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

intervals = blast.windowmasker_mask_sequences("genome.fa", "counts.obin", "mask.int", verbose=False)
soft = blast.windowmasker_mask_sequences("genome.fa", "counts.txt", "soft.fa", output_format="fasta", verbose=False)
lines = open("mask.int").read().splitlines()
lowered = sum(char.islower() for char in open("soft.fa").read().split("\n", 1)[1])
print(intervals["status"], soft["status"], lines[0], lowered > 0)
```

- **示例输出：**

```text
completed completed >chr1 True
```

## 8) `dnakit.alignment.convert2blastmask_from_fasta()`

- **作用：** 把 FASTA 中的小写屏蔽转换为 `makeblastdb -mask_data` 可用的屏蔽文件。
- **处理方法：** 运行 `convert2blastmask`，传入屏蔽算法名称和对其选项的自由文本说明。
- **API：** `dnakit.alignment.convert2blastmask_from_fasta(input_path, output, *, masking_algorithm, masking_options, output_format='maskinfo_asn1_text', parse_seqids=False, blast_path=None)`.
- **输入：** 软屏蔽 FASTA、新的输出路径、`masking_algorithm`（`dust`、`seg`、`windowmasker`、`repeat`、`other`）、`masking_options`，以及可选的 `output_format`、`parse_seqids`。
- **输出：** 通用结果键；屏蔽文件即输出。
- **示例代码：**

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.convert2blastmask_from_fasta("soft.fa", "mask.asn", masking_algorithm="windowmasker", masking_options="-mk_counts", verbose=False)
print(result["status"], open("mask.asn").read().split(" ::=")[0])
```

- **示例输出：**

```text
completed Blast-db-mask-info
```

## 9) `dnakit.alignment.blast_vdb_search()`

- **作用：** 用核酸查询（`blastn_vdb`）或蛋白查询（`tblastn_vdb`）搜索 SRA 运行号。
- **处理方法：** 运行 `blastn_vdb` 或 `tblastn_vdb`，把运行号作为数据库。程序通过网络从 NCBI 读取这些运行数据。
- **API：** `dnakit.alignment.blast_vdb_search(query, accessions, output, *, program='blastn_vdb', outfmt=None, evalue=None, max_target_seqs=None, num_threads=None, task=None, genetic_code=None, sra_mode=None, include_filtered_reads=False, blast_path=None)`.
- **输入：** 查询 FASTA、SRA 运行号列表、新的输出路径、`program`，以及可选的 `outfmt`、`evalue`、`max_target_seqs`、`num_threads`、`task`（blastn_vdb）、`genetic_code`（tblastn_vdb）、`sra_mode`、`include_filtered_reads`。
- **输出：** 通用结果键；表格格式时还有 `rows` / `row_count`。
- **示例调用（此处不运行）：** 该调用会从 NCBI 下载 SRA 数据，因此不在文档测试中运行。

```python
from dnakit.alignment import blast_p14_cli_functions as blast

result = blast.blast_vdb_search("read.fa", ["SRR000001"], "vdb_hits.tsv", outfmt=6, max_target_seqs=5)
print(result["status"], result["row_count"])
```

## 10) `dnakit.alignment.blast_vdb_cmd_fetch_reads()`

- **作用：** 按 SRA 运行号取出序列记录、运行信息或参考序列。
- **处理方法：** 运行 `blast_vdb_cmd`。`view='sequences'` 写出 `entries`（或 `['all']`）或 `entry_batch` 指定的 reads；`'info'` 输出每个运行的序列数和碱基数；`'references'` 导出已比对运行的参考序列；`'paths'`、`'paths_all'` 列出运行文件。运行数据通过网络从 NCBI 读取。
- **API：** `dnakit.alignment.blast_vdb_cmd_fetch_reads(accessions, output, *, view='sequences', entries=None, entry_batch=None, outfmt=None, line_length=None, include_filtered_reads=False, blast_path=None)`.
- **输入：** SRA 运行号列表、新的输出路径、`view`；`sequences` 时需 `entries` 与 `entry_batch` 二选一，另有可选的 `outfmt`、`line_length`、`include_filtered_reads`。
- **输出：** 通用结果键；记录写在输出文件中。
- **示例调用（此处不运行）：** 该调用需要联网访问 NCBI，因此不在文档测试中运行。注意 `blast_vdb_cmd` 在 NCBI 的 BLAST+ 压缩包中，但不在 bioconda 的 `blast` 包中。

```python
from dnakit.alignment import blast_p14_cli_functions as blast

info = blast.blast_vdb_cmd_fetch_reads(["SRR000001"], "run_info.txt", view="info")
print(info["status"])
print(open("run_info.txt").read())
```

## 11) `dnakit.alignment.bamtools_convert_alignments()`

- **作用：** 把 BAM 比对转换为 JSON、YAML 或其他文本格式。
- **处理方法：** 运行 `bamtools convert`，`-format` 取 `json`（每条比对一个 JSON 对象：名称、flag、参考序列、位置、比对质量、CIGAR、碱基、质量值、标签）、`yaml`、`bed`、`fasta`、`fastq`、`pileup` 或 `sam`。
- **API：** `dnakit.alignment.bamtools_convert_alignments(bams, output, *, output_format='json', region=None, no_header=False, reference_fasta=None, print_mapping_quality=False, bamtools_path=None)`.
- **输入：** BAM 文件列表（会合并）、新的输出路径、`output_format`，以及可选的 `region`、`no_header`（sam）、`reference_fasta` 和 `print_mapping_quality`（pileup）。
- **输出：** 通用结果键；转换后的记录写在输出文件中。
- **示例代码：**

```python
import json

from dnakit.alignment import bamtools_p14_cli_functions as bt

bt.bamtools_convert_alignments(["reads.bam"], "chr2.json", output_format="json", region="chr2", verbose=False)
records = [json.loads(line) for line in open("chr2.json") if line.strip()]
bt.bamtools_convert_alignments(["reads.bam"], "chr2.yaml", output_format="yaml", region="chr2", verbose=False)
print(len(records), records[0]["reference"], sorted(records[0])[:4], open("chr2.yaml").read().count("RefName: chr2"))
```

- **示例输出：**

```text
20 chr2 ['alignmentFlag', 'cigar', 'filename', 'mapQuality'] 20
```

## 12) `dnakit.alignment.bamtools_count_alignments()`

- **作用：** 统计一个或多个 BAM 文件中的比对数量，可限定区域。
- **处理方法：** 运行 `bamtools count`；文件先合并再计数，返回程序打印的数字。
- **API：** `dnakit.alignment.bamtools_count_alignments(bams, *, region=None, bamtools_path=None)`.
- **输入：** BAM 文件列表和可选的 `region`。
- **输出：** 通用结果键，加上 `count`。
- **示例代码：**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

print(bt.bamtools_count_alignments(["reads.bam"], verbose=False)["count"], bt.bamtools_count_alignments(["reads.bam"], region="chr1", verbose=False)["count"])
```

- **示例输出：**

```text
60 20
```

## 13) `dnakit.alignment.bamtools_filter_alignments()`

- **作用：** 保留满足全部给定属性条件的 BAM 比对。
- **处理方法：** 运行 `bamtools filter`，条件由参数构成：精确 flag、插入长度/读长/比对质量的比较式（`50`、`>=50`、`<100`、`!0`）、名称（可含 `*`）、碱基模式、标签值（`NM:<4`）以及 12 个 `is_*` 标志过滤，或使用过滤脚本。
- **API：** `dnakit.alignment.bamtools_filter_alignments(bams, output, *, region=None, script=None, alignment_flag=None, insert_size=None, length=None, map_quality=None, name=None, query_bases=None, tag=None, is_duplicate=None, is_failed_qc=None, is_first_mate=None, is_mapped=None, is_mate_mapped=None, is_mate_reverse_strand=None, is_paired=None, is_primary_alignment=None, is_proper_pair=None, is_reverse_strand=None, is_second_mate=None, is_singleton=None, bamtools_path=None)`.
- **输入：** BAM 文件列表、新的输出路径、至少一个条件或 `script`，以及可选的 `region`。
- **输出：** 通用结果键；保留的比对写在输出 BAM 中。
- **示例代码：**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_filter_alignments(["reads.bam"], "filtered.bam", map_quality=">=60", is_reverse_strand=False, is_mapped=True, verbose=False)
print(result["status"], bt.bamtools_count_alignments(["filtered.bam"], verbose=False)["count"])
```

- **示例输出：**

```text
completed 19
```

## 14) `dnakit.alignment.bamtools_index_bam_bti()`

- **作用：** 为按坐标排序的 BAM 文件建立 BamTools 专用索引（`.bti`）。
- **处理方法：** 运行 `bamtools index -bti`，不生成标准 BAI 索引。索引写在 BAM 旁边，且该文件必须尚不存在。
- **API：** `dnakit.alignment.bamtools_index_bam_bti(bam, *, bamtools_path=None)`.
- **输入：** 已排序的 BAM 文件。
- **输出：** 通用结果键；`outputs` 为 `.bti` 路径。
- **示例代码：**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_index_bam_bti("reads.bam", verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **示例输出：**

```text
completed ['reads.bam.bti']
```

## 15) `dnakit.alignment.bamtools_sample_alignments()`

- **作用：** 从一个或多个 BAM 文件中随机抽取一部分比对写出。
- **处理方法：** 运行 `bamtools random`，带 `-n` 和 `-seed`。程序不做重复检查，同一条比对可能被抽到两次。
- **API：** `dnakit.alignment.bamtools_sample_alignments(bams, output, *, count=None, seed=None, region=None, bamtools_path=None)`.
- **输入：** BAM 文件列表、新的输出路径，以及可选的 `count`（BamTools 默认 10000）、`seed`、`region`。
- **输出：** 通用结果键；抽样结果写在输出 BAM 中。
- **示例代码：**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_sample_alignments(["reads.bam"], "sample.bam", count=15, seed=1, verbose=False)
print(result["status"], bt.bamtools_count_alignments(["sample.bam"], verbose=False)["count"])
```

- **示例输出：**

```text
completed 15
```

## 16) `dnakit.alignment.bamtools_resolve_pairs()`

- **作用：** 根据片段长度模型标记双端比对的 proper-pair 标志。
- **处理方法：** 运行 `bamtools resolve`。`mode='twoPass'` 一次完成统计学习与标记；`'makeStats'` 只写统计文件；`'markPairs'` 用已有统计文件做标记（BamTools 还需要同目录的 `<stats>.uniq_names.txt`）。
- **API：** `dnakit.alignment.bamtools_resolve_pairs(bam, stats, *, output=None, mode='twoPass', min_mapping_quality=None, confidence_interval=None, unused_model_threshold=None, force=False, bamtools_path=None)`.
- **输入：** BAM 文件、`stats` 路径（twoPass/makeStats 为新文件，markPairs 为已有文件）、`output` BAM（makeStats 不需要）、`mode`，以及可选的 `min_mapping_quality`、`confidence_interval`、`unused_model_threshold`、`force`。
- **输出：** 通用结果键，加上 `mode`。
- **示例代码：**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

stats = bt.bamtools_resolve_pairs("reads.bam", "fragments.txt", mode="makeStats", verbose=False)
marked = bt.bamtools_resolve_pairs("reads.bam", "fragments.txt", output="resolved.bam", mode="markPairs", verbose=False)
print(stats["status"], marked["status"], marked["mode"])
```

- **示例输出：**

```text
completed completed markPairs
```

## 17) `dnakit.alignment.bamtools_revert_alignments()`

- **作用：** 去除 BAM 比对的重复标记，并恢复原始碱基质量值（OQ 标签）。
- **处理方法：** 运行 `bamtools revert`；`keep_duplicate_marks` 和 `keep_qualities` 分别关闭这两项操作。
- **API：** `dnakit.alignment.bamtools_revert_alignments(bams, output, *, keep_duplicate_marks=False, keep_qualities=False, bamtools_path=None)`.
- **输入：** BAM 文件列表、新的输出路径，以及可选的 `keep_duplicate_marks`、`keep_qualities`。
- **输出：** 通用结果键；还原后的比对写在输出 BAM 中。
- **示例代码：**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

result = bt.bamtools_revert_alignments(["reads.bam"], "reverted.bam", verbose=False)
print(result["status"], bt.bamtools_count_alignments(["reverted.bam"], verbose=False)["count"])
```

- **示例输出：**

```text
completed 60
```

## 18) `dnakit.alignment.bamtools_split_alignments()`

- **作用：** 按某个属性的取值把 BAM 文件拆分为多个 BAM 文件。
- **处理方法：** 运行 `bamtools split`，用 `-mapped`（`.MAPPED.bam`、`.UNMAPPED.bam`）、`-paired`（`.PAIRED_END.bam`、`.SINGLE_END.bam`）、`-reference`（`.REF_<名称>.bam`）或 `-tag <TAG>`（`.TAG_<标签>_<取值>.bam`）。只有实际出现的取值才会生成文件。
- **API：** `dnakit.alignment.bamtools_split_alignments(bam, output_stub, *, by='mapped', tag=None, reference_prefix=None, tag_prefix=None, tag_list_delimiter=None, bamtools_path=None)`.
- **输入：** BAM 文件、输出前缀 `output_stub`、`by`（`mapped`、`paired`、`reference`、`tag`）、按标签拆分时的 `tag`，以及可选的 `reference_prefix`、`tag_prefix`、`tag_list_delimiter`。
- **输出：** 通用结果键；`outputs` 列出写出的文件。
- **示例代码：**

```python
from dnakit.alignment import bamtools_p14_cli_functions as bt

for by, tag in (("mapped", None), ("reference", None), ("tag", "XX")):
    result = bt.bamtools_split_alignments("reads.bam", f"split_{by}", by=by, tag=tag, verbose=False)
    print(by, [path.rsplit(".", 2)[-2] for path in result["outputs"]])
```

- **示例输出：**

```text
mapped ['MAPPED', 'UNMAPPED']
reference ['REF_chr1', 'REF_chr2', 'REF_unmapped']
tag ['TAG_XX_a', 'TAG_XX_b']
```

## 19) `dnakit.alignment.hisat2_extract_splice_sites()`

- **作用：** 从 GTF 基因注释中提取剪接位点。
- **处理方法：** 运行 HISAT2 辅助脚本 `hisat2_extract_splice_sites.py`。输出每行对应一个内含子：染色体、上游外显子最后一个碱基、下游外显子第一个碱基（0 起始）和链。
- **API：** `dnakit.alignment.hisat2_extract_splice_sites(annotation_gtf, output, *, hisat2_path=None)`.
- **输入：** GTF 文件和新的输出路径。
- **输出：** 通用结果键，加上 `rows` / `row_count`。
- **示例代码：**

```python
import random

from dnakit.alignment import hisat2_p14_cli_functions as hs

rng = random.Random(8)
sequence = "".join(rng.choice("ACGT") for _ in range(3000))
open("ref.fa", "w").write(">chr1\n" + "\n".join(sequence[i:i + 60] for i in range(0, 3000, 60)) + "\n")
exons = [(201, 400), (601, 800), (1001, 1200)]
open("genes.gtf", "w").write("".join(f'chr1\tdemo\texon\t{a}\t{b}\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n' for a, b in exons))
result = hs.hisat2_extract_splice_sites("genes.gtf", "ss.txt", verbose=False)
print(result["status"], result["rows"])
```

- **示例输出：**

```text
completed [['chr1', '399', '600', '+'], ['chr1', '799', '1000', '+']]
```

## 20) `dnakit.alignment.hisat2_extract_exons()`

- **作用：** 从 GTF 基因注释中提取外显子。
- **处理方法：** 运行 `hisat2_extract_exons.py`。输出每行是一个合并后的外显子：染色体、起点、终点（0 起始）和链，可作为 `hisat2-build --exon` 的输入。
- **API：** `dnakit.alignment.hisat2_extract_exons(annotation_gtf, output, *, hisat2_path=None)`.
- **输入：** GTF 文件和新的输出路径。
- **输出：** 通用结果键，加上 `rows` / `row_count`。
- **示例代码：**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

result = hs.hisat2_extract_exons("genes.gtf", "exons.txt", verbose=False)
print(result["status"], result["row_count"], result["rows"][0])
```

- **示例输出：**

```text
completed 3 ['chr1', '200', '399', '+']
```

## 21) `dnakit.alignment.hisat2_inspect_index()`

- **作用：** 输出 HISAT2 索引中保存的 SNP、剪接位点或外显子。
- **处理方法：** 运行 `hisat2-inspect`，带 `--snp`、`--ss`、`--ss-all` 或 `--exon`（`view` 取 `snps`、`splice_sites`、`splice_sites_all`、`exons`）；`summary`、`names`、`sequences` 分别给出索引概要、参考序列名和被索引的 FASTA。
- **API：** `dnakit.alignment.hisat2_inspect_index(index_prefix, output, *, view='snps', hisat2_path=None)`.
- **输入：** `hisat2-build` 生成的索引前缀、新的输出路径和 `view`。
- **输出：** 通用结果键，加上 `view`；四种注释视图还有 `rows` / `row_count`。
- **示例代码：**

```python
import subprocess

from dnakit.alignment import hisat2_p14_cli_functions as hs

subprocess.run(["hisat2-build", "--ss", "ss.txt", "--exon", "exons.txt", "ref.fa", "ref_idx"], check=True, capture_output=True)
for view in ("splice_sites", "exons"):
    result = hs.hisat2_inspect_index("ref_idx", f"inspect_{view}.txt", view=view, verbose=False)
    print(view, result["status"], result["row_count"])
```

- **示例输出：**

```text
splice_sites completed 2
exons completed 3
```

## 22) `dnakit.alignment.hisat2_align_reads()`

- **作用：** 用支持剪接比对的 HISAT2 把 reads 比对到已建索引的基因组。
- **处理方法：** 运行 `hisat2 -x <索引>`，输入单端（`-U`）或双端（`-1`/`-2`）reads，输出 SAM。`no_spliced_alignment` 用于 DNA reads；`known_splice_sites`、`novel_splice_sites_out`、`rna_strandness`、`downstream_transcriptome_assembly`（`--dta`）等选项对应含义相同的 HISAT2 参数。
- **API：** `dnakit.alignment.hisat2_align_reads(index_prefix, output_sam, *, unpaired_reads=None, reads_1=None, reads_2=None, read_format='fastq', preset=None, threads=None, no_spliced_alignment=False, known_splice_sites=None, novel_splice_sites_out=None, novel_splice_sites_in=None, min_intron_length=None, max_intron_length=None, rna_strandness=None, downstream_transcriptome_assembly=False, dta_cufflinks=False, transcriptome_mapping_only=False, max_alignments=None, min_fragment_length=None, max_fragment_length=None, no_mixed=False, no_discordant=False, trim_5_prime=None, trim_3_prime=None, skip_reads=None, stop_after_reads=None, read_group_id=None, no_softclip=False, no_repeat_index=False, summary_file=None, seed=None, hisat2_path=None)`.
- **输入：** 索引前缀、新的 SAM 路径、`unpaired_reads` 和/或 `reads_1` + `reads_2`，以及可选的 `read_format`、`preset`、`threads`、`known_splice_sites`、`summary_file`、`max_alignments` 和签名中列出的其他选项。
- **输出：** 通用结果键；写出 SAM 文件（以及可选的汇总和剪接位点文件）。
- **示例代码：**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

sequence = "".join(line.strip() for line in open("ref.fa").read().split("\n")[1:])
open("reads.fq", "w").write(
    f"@spliced\n{sequence[350:400] + sequence[600:650]}\n+\n{'I' * 100}\n@plain\n{sequence[1500:1600]}\n+\n{'I' * 100}\n"
)
result = hs.hisat2_align_reads("ref_idx", "aligned.sam", unpaired_reads=["reads.fq"], known_splice_sites="ss.txt", verbose=False)
rows = {line.split("\t")[0]: line.split("\t") for line in open("aligned.sam") if not line.startswith("@")}
print(result["status"], rows["plain"][3], rows["plain"][5], rows["spliced"][5])
```

- **示例输出：**

```text
completed 1501 100M 50M200N50M
```

## 23) `dnakit.alignment.hisat2_read_statistics()`

- **作用：** 统计 reads 数量，并报告最短、最长和平均读长。
- **处理方法：** 对前 `read_count` 条 reads 运行 `hisat2_read_statistics.py`（脚本默认 10000）。
- **API：** `dnakit.alignment.hisat2_read_statistics(reads, *, read_count=None, hisat2_path=None)`.
- **输入：** reads 文件和可选的 `read_count`。
- **输出：** 通用结果键，加上 `reads_examined`、`minimum_length`、`maximum_length`、`average_length`、`read_lengths`。
- **示例代码：**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

result = hs.hisat2_read_statistics("reads.fq", verbose=False)
print(result["reads_examined"], result["minimum_length"], result["maximum_length"], result["average_length"])
```

- **示例输出：**

```text
2 100 100 100
```

## 24) `dnakit.alignment.hisat2_extract_snps_haplotypes_vcf()`

- **作用：** 从 VCF 文件提取 SNP 和单倍型，供 `hisat2-build` 使用。
- **处理方法：** 运行 `hisat2_extract_snps_haplotypes_VCF.py`，写出 `<base>.snp` 和 `<base>.haplotype`（`extra_files` 时还有额外文件）。
- **API：** `dnakit.alignment.hisat2_extract_snps_haplotypes_vcf(genome_fasta, vcf_files, output_base, *, reference_type=None, inter_gap=None, intra_gap=None, non_rs_ids=False, genotype_vcf=None, genotype_gene_list=None, extra_files=False, hisat2_path=None)`.
- **输入：** 基因组 FASTA、VCF 文件列表、输出基础路径，以及可选的 `reference_type`、`inter_gap`、`intra_gap`、`non_rs_ids`、`genotype_vcf`、`genotype_gene_list`、`extra_files`。
- **输出：** 通用结果键；`outputs` 列出写出的文件。
- **示例代码：**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

sequence = "".join(line.strip() for line in open("ref.fa").read().split("\n")[1:])
flip = {"A": "G", "C": "T", "G": "A", "T": "C"}
rows = "".join(f"chr1\t{pos}\trs{pos}\t{sequence[pos - 1]}\t{flip[sequence[pos - 1]]}\t.\t.\t.\n" for pos in (200, 500, 900))
open("variants.vcf", "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + rows)
result = hs.hisat2_extract_snps_haplotypes_vcf("ref.fa", ["variants.vcf"], "vcf_snps", verbose=False)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
print(open("vcf_snps.snp").read().splitlines()[0].split("\t")[1:4])
```

- **示例输出：**

```text
completed ['vcf_snps.haplotype', 'vcf_snps.snp']
['single', 'chr1', '199']
```

## 25) `dnakit.alignment.hisat2_extract_snps_haplotypes_ucsc()`

- **作用：** 从 UCSC dbSNP 表提取 SNP 和单倍型，供 `hisat2-build` 使用。
- **处理方法：** 对 UCSC SNP 表（如 `snp144Common.txt`，纯文本或 gzip）运行 `hisat2_extract_snps_haplotypes_UCSC.py`，写出 `<base>.snp` 和 `<base>.haplotype`。
- **API：** `dnakit.alignment.hisat2_extract_snps_haplotypes_ucsc(genome_fasta, snp_table, output_base, *, inter_gap=None, intra_gap=None, hisat2_path=None)`.
- **输入：** 基因组 FASTA、UCSC SNP 表、输出基础路径，以及可选的 `inter_gap`、`intra_gap`。
- **输出：** 通用结果键；`outputs` 列出写出的文件。
- **示例代码：**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

sequence = "".join(line.strip() for line in open("ref.fa").read().split("\n")[1:])
flip = {"A": "G", "C": "T", "G": "A", "T": "C"}
lines = []
for index, pos in enumerate((100, 300, 900)):
    ref = sequence[pos]
    fields = ["585", "chr1", str(pos), str(pos + 1), f"rs{1000 + index}", "0", "+", ref, ref, f"{ref}/{flip[ref]}", "genomic", "single",
              "unknown", "0", "0", "unknown", "exact", "1", "", "1", "", "2", f"{ref},{flip[ref]},", "10,2,", "0.8,0.2,", "0"]
    lines.append("\t".join(fields))
open("snp_table.txt", "w").write("\n".join(lines) + "\n")
result = hs.hisat2_extract_snps_haplotypes_ucsc("ref.fa", "snp_table.txt", "ucsc_snps", verbose=False)
print(result["status"], len(open("ucsc_snps.snp").read().splitlines()))
```

- **示例输出：**

```text
completed 3
```

## 26) `dnakit.alignment.hisat2_simulate_reads()`

- **作用：** 根据基因组、注释和 SNP 模拟 RNA-seq reads（默认）或 DNA-seq reads。
- **处理方法：** 运行 `hisat2_simulate_reads.py`；`dna=True` 时加 `--dna`。该脚本只保留名为 1-22、X、Y 的染色体，其 RNA-seq 模式需要 Python 3.10 或更早版本（调用了 3.11 已移除的 `random.shuffle(..., random=...)`）；可用 `python_path` 指定旧版解释器。
- **API：** `dnakit.alignment.hisat2_simulate_reads(genome_fasta, annotation_gtf, snp_file, output_base, *, dna=False, single_end=False, read_length=None, fragment_length=None, num_fragments=None, expression_profile=None, repeat_info=None, error_rate=None, max_mismatch=None, random_seed=None, snp_probability=None, python_path=None, hisat2_path=None)`.
- **输入：** 基因组 FASTA、GTF、SNP 文件、输出基础路径、`dna`、`single_end`，以及可选的 `read_length`、`fragment_length`、`num_fragments`、`expression_profile`、`repeat_info`、`error_rate`、`max_mismatch`、`random_seed`、`snp_probability`、`python_path`。
- **输出：** 通用结果键，加上 `read_type`；`outputs` 列出 reads 与真值文件。
- **示例代码：**

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

open("ref_num.fa", "w").write(open("ref.fa").read().replace(">chr1", ">1"))
open("genes_num.gtf", "w").write(open("genes.gtf").read().replace("chr1", "1"))
open("snps_num.snp", "w").write(open("vcf_snps.snp").read().replace("chr1", "1"))
result = hs.hisat2_simulate_reads("ref_num.fa", "genes_num.gtf", "snps_num.snp", "sim", dna=True, single_end=True,
                                  num_fragments=50, read_length=40, random_seed=1, verbose=False)
print(result["status"], result["read_type"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```

- **示例输出：**

```text
completed dna ['sim.sam', 'sim_1.fa']
```

## 27) `dnakit.alignment.hisat2_build_repeat_index()`

- **作用：** 在参考基因组中寻找重复序列，并写出 HISAT2 重复索引所需的文件。
- **处理方法：** 运行 `hisat2-repeat`，带种子、重复长度、重复拷贝数和编辑距离选项。该程序属于 HISAT2 源码构建，不在 bioconda 软件包中。
- **API：** `dnakit.alignment.hisat2_build_repeat_index(reference, output_prefix, *, threads=None, seed_length=None, seed_count=None, min_repeat_length=None, max_repeat_length=None, repeat_count=None, repeat_edit_distance=None, forward_only=False, hisat2_path=None)`.
- **输入：** 参考 FASTA、新的输出前缀，以及可选的 `threads`、`seed_length`、`seed_count`、`min_repeat_length`、`max_repeat_length`、`repeat_count`、`repeat_edit_distance`、`forward_only`。
- **输出：** 通用结果键；`outputs` 列出 `.rep.*` 文件。
- **示例调用（此处不运行）：** `hisat2-repeat` 需从 HISAT2 源码编译（`make hisat2-repeat`），因此该调用不在文档测试中运行；程序在 PATH 中时单元测试会运行它。

```python
from dnakit.alignment import hisat2_p14_cli_functions as hs

result = hs.hisat2_build_repeat_index("genome.fa", "genome_repeats", seed_length=30)
print(result["status"], [path.rsplit("/", 1)[-1] for path in result["outputs"]])
```
