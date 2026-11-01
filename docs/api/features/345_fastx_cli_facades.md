# FASTA/FASTQ 命令行适配器：seqmagick、SeqKit、AfterQC 与 Trimmomatic

这一组函数通过固定的参数白名单、不经过 shell 调用四个由使用者自行安装的程序：[seqmagick](https://github.com/fhcrc/seqmagick)（按位置切片残基、比对清理、记录筛选、质量过滤、文件汇总）、[SeqKit](https://github.com/shenwei356/seqkit)（FASTA 索引、修复损坏文件、实时汇集、共有记录、序列名替换、BAM 工具箱）、[AfterQC](https://github.com/OpenGene/AfterQC)（重叠分析、气泡去除、条形码）和 [Trimmomatic](https://github.com/usadellab/Trimmomatic)（修剪与过滤步骤、配对拆分）。这些程序不随包提供：缺少程序时报错并给出安装命令。输出文件必须尚不存在，不会覆盖任何文件。每个结果都是字典，含 `status`、`backend`、`version`、`command`、`outputs`、解析后的表和 `log_tail`。

测试版本：seqmagick 0.8.6、SeqKit 2.14.0、Trimmomatic 0.39。AfterQC 0.9.7 是 Python 2 程序，测试环境中没有，因此它的三个函数只测试了参数处理。示例运行 seqmagick、SeqKit 和 Trimmomatic。

## 示例数据

所有输入都是很小的文本文件；示例把结果写入临时文件夹。

**示例文件：** `aln.fa`

```text
>s1
AC--GTACGT
>s2
AC--GTTCGT
>s3
ACGTGTACGA
>s4
----------
```

**示例文件：** `x.fa`

```text
>a first
ACGTACGTAC
>b second
GGGGCCCCAA
>c third
TTTTTTTTTT
```

**示例文件：** `other.fa`

```text
>a other
ACGTACGTAC
>x only
CCCCCCCCCC
>c third
AAAAAAAAAA
```

**示例文件：** `reads.fq`

```text
@r1
CAGATGACGTACGTACGT
+
IIIIIIIIIIIIIIIIII
@r2
TTTTACGTACGTACGTAC
+
IIIIIIIIIIIIIIIIII
@r3
CAGATGACNTACGTACGT
+
IIIIIIIIIIIIIIIIII
```

**示例文件：** `broken.fq`

```text
@r1
ACGT
+
IIII
@bad
AC
+
IIII
@r3
ACGTN
+
IIIII
nonsense
```

**示例文件：** `r1.fq`

```text
@p1/1
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@p2/1
TTGCATTGCATTGCATTGCATTGCATTGCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**示例文件：** `r2.fq`

```text
@p1/2
GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@p2/2
ACACAACACAACACAACACAACACAACACA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

**示例文件：** `m1.fq`

```text
@m1/1
ACGTAC
+
IIIIII
@m2/1
ACGTAC
+
IIIIII
@m3/1
ACGTAC
+
IIIIII
```

**示例文件：** `m2.fq`

```text
@m1/2
GGGGGG
+
IIIIII
@m3/2
TTTTTT
+
IIIIII
```

## 1) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_slice_residues()`

- **作用：** 按位置区间保留、删除或遮盖每条序列的残基。
- **处理方法：** 运行 `seqmagick convert --cut/--drop/--mask`。区间为 1 起始、含两端的 `start:end`，多个区间用逗号分隔；`relative_to` 按某条含空位序列的非空位残基编号。
- **API：** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_slice_residues(source, destination, *, cut=None, drop=None, mask=None, relative_to=None, input_format=None, output_format=None, alphabet=None, line_wrap=None, seqmagick_path=None)`。
- **输入：** 序列文件、新的目标路径，以及 `cut`、`drop`、`mask` 中至少一个。
- **输出：** `output`、`output_record_count`、命令与版本。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_slice_residues

result = seqmagick_slice_residues("aln.fa", str(out / "cut.fa"), cut="3:6")
print(result["output_record_count"], result["version"])
print((out / "cut.fa").read_text().splitlines()[1])
```

- **示例输出：**

```text
4 0.8.6
--GT
```

## 2) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_clean_alignment()`

- **作用：** 统一空位字符、去掉以空位为主的列并删除全空位序列。
- **处理方法：** 运行 `seqmagick convert --dash-gap`、`--squeeze` / `--squeeze-threshold` 与 `--prune-empty`。
- **API：** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_clean_alignment(source, destination, *, dash_gap=False, squeeze=False, squeeze_threshold=None, prune_empty=False, input_format=None, output_format=None, alphabet=None, line_wrap=None, seqmagick_path=None)`。
- **输入：** 比对文件、新的目标路径和至少一个选项。
- **输出：** `output`、`output_record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_clean_alignment

result = seqmagick_clean_alignment("aln.fa", str(out / "clean.fa"), squeeze_threshold=0.75, prune_empty=True)
print(result["output_record_count"])
print((out / "clean.fa").read_text().splitlines()[1])
```

- **示例输出：**

```text
3
ACGTACGT
```

## 3) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_select_sequences()`

- **作用：** 按 ID 列表或序列内容的正则表达式保留或删除序列。
- **处理方法：** 运行带 `--include-from-file`、`--exclude-from-file`、`--seq-pattern-include` 或 `--seq-pattern-exclude` 的 `seqmagick convert`。
- **API：** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_select_sequences(source, destination, *, include_ids_file=None, exclude_ids_file=None, seq_pattern_include=None, seq_pattern_exclude=None, input_format=None, output_format=None, alphabet=None, line_wrap=None, seqmagick_path=None)`。
- **输入：** 序列文件、新的目标路径和至少一个筛选条件。
- **输出：** `output`、`output_record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_select_sequences

ids = out / "ids.txt"
ids.write_text("s1\ns3\n")
result = seqmagick_select_sequences("aln.fa", str(out / "kept.fa"), include_ids_file=str(ids))
print(result["output_record_count"])
result = seqmagick_select_sequences("aln.fa", str(out / "pattern.fa"), seq_pattern_include="GTAC")
print(result["output_record_count"])
```

- **示例输出：**

```text
2
2
```

## 4) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_quality_filter()`

- **作用：** 按质量、长度、模糊碱基和 5' 端 IUPAC 引物过滤 reads。
- **处理方法：** 运行 `seqmagick quality-filter`，默认值与原程序一致（平均质量 25，长度 200 至 1000）。`ambiguous_action` 在 `N` 处截断或丢弃 reads。引物通过 seqmagick 的条形码机制检查，因此函数提供一行空条形码的表；通过的 reads 会去掉引物。
- **API：** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_quality_filter(sequence_file, output_file, *, input_qual=None, min_mean_quality=25.0, min_length=200, max_length=1000, ambiguous_action=None, max_ambiguous=None, pct_ambiguous=None, primer=None, barcode_file=None, barcode_header=False, sample_map_out=None, quality_window=0, quality_window_mean_qual=None, quality_window_prop=1.0, details_out=None, seqmagick_path=None)`。
- **输入：** FASTQ（或 FASTA 加 `input_qual`）文件、新的输出文件、阈值，以及可选的 `primer`、`barcode_file`。
- **输出：** `output`、`output_record_count`、`report`（每个过滤器一行，含通过与失败数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_quality_filter

result = seqmagick_quality_filter("reads.fq", str(out / "q.fq"), min_length=5, max_length=100, ambiguous_action="truncate")
print(result["output_record_count"])
for row in result["report"]:
    print(row["name"], row["failed"])
primed = seqmagick_quality_filter("reads.fq", str(out / "p.fq"), min_length=5, max_length=100, primer="CAGATG")
print(primed["output_record_count"], (out / "p.fq").read_text().splitlines()[1])
```

- **示例输出：**

```text
3
Quality Score [min_mean: 25.0] 0
Maximum Length [100] 0
Minimum Length [5] 0
Ambiguous Base [truncate] 0
2 ACGTACGTACGT
```

## 5) `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_info()`

- **作用：** 汇总序列文件：记录数、长度范围、平均长度以及是否为比对。
- **处理方法：** 运行 `seqmagick info`；所有序列等长时该文件视为比对。
- **API：** `dnakit.fasta_fastq.seqmagick_p15_cli_functions.seqmagick_info(sequence_files, output_file, *, input_format=None, table_format='tab', threads=1, seqmagick_path=None)`。
- **输入：** 文件列表、新的输出文件和 `table_format`（`tab`、`csv`、`align`）。
- **输出：** `files`：每个文件一行，含 `alignment`、`min_len`、`max_len`、`avg_len`、`num_seqs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqmagick_p15_cli_functions import seqmagick_info

result = seqmagick_info(["aln.fa", "reads.fq"], str(out / "info.tsv"))
for row in result["files"]:
    print(Path(row["name"]).name, row["alignment"], row["min_len"], row["max_len"], row["num_seqs"])
```

- **示例输出：**

```text
aln.fa TRUE 10 10 4
reads.fq TRUE 18 18 3
```

## 6) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_faidx()`

- **作用：** 建立 FASTA 索引（`.fai`），并可提取区域。
- **处理方法：** 运行 `seqkit faidx`。区域为 1 起始，可用负数位置；`full_head` 生成 `.seqkit.fai` 并输出完整标题行。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_faidx(fasta_path, *, regions=None, output_path=None, full_head=False, use_regexp=False, ignore_case=False, update_index=True, seqkit_path=None)`。
- **输入：** FASTA 文件；可选 `regions` 并配新的 `output_path`。
- **输出：** `index_file`、`sequence_count`、`index` 行（名称、长度、偏移、每行碱基数、行宽）、`extracted_records`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import shutil
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_faidx

fasta = shutil.copy("x.fa", out / "x.fa")
result = seqkit_faidx(str(fasta))
print(result["sequence_count"], result["index"][1])
region = seqkit_faidx(str(fasta), regions=["a:2-4", "b"], output_path=str(out / "regions.fa"))
print(region["extracted_records"], (out / "regions.fa").read_text().splitlines()[1])
```

- **示例输出：**

```text
3 {'name': 'b', 'length': 10, 'offset': 30, 'line_bases': 10, 'line_width': 11}
2 CGT
```

## 7) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_sanitize_fastx()`

- **作用：** 跳过格式错误的记录，修复损坏的单行 FASTQ/FASTA 文件。
- **处理方法：** 运行 `seqkit sana`，出错后继续解析，丢弃序列与质量长度不一致的记录。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_sanitize_fastx(input_path, output_path, *, format='fastq', input_format=None, output_format=None, allow_gaps=False, quality_ascii_base=33, seqkit_path=None)`。
- **输入：** 损坏的文件和新的输出路径；`format`、`input_format`、`output_format`。
- **输出：** `passed_records`、`discarded_lines`、`output_record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_sanitize_fastx

result = seqkit_sanitize_fastx("broken.fq", str(out / "clean.fq"))
print(result["passed_records"], result["discarded_lines"], result["output_record_count"])
```

- **示例输出：**

```text
2 5 2
```

## 8) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_stream_directory()`

- **作用：** 实时汇集并输出目录中陆续出现的 FASTX 文件。
- **处理方法：** 运行 `seqkit scat`。`find_only` 只合并已有文件并结束；否则必须给出 `time_limit`，保证调用会结束。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_stream_directory(directory, output_path, *, format='fastq', find_only=False, time_limit=None, gz_only=False, file_pattern=None, minimum_size_increase_kb=5, quality_ascii_base=33, seqkit_path=None)`。
- **输入：** 目录、新的输出路径、`format`，以及 `find_only` 或 `time_limit`。
- **输出：** `passed_records`、`output_record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import shutil
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_stream_directory

run = out / "run"
run.mkdir()
shutil.copy("r1.fq", run / "one.fq")
shutil.copy("r2.fq", run / "two.fastq")
result = seqkit_stream_directory(str(run), str(out / "all.fq"), find_only=True)
print(result["passed_records"], result["output_record_count"])
```

- **示例输出：**

```text
4 4
```

## 9) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_common_records()`

- **作用：** 两个或多个文件共有的记录，按 ID、完整名称或序列匹配。
- **处理方法：** 运行 `seqkit common`。按序列匹配时默认比较正反两条链，除非设置 `only_positive_strand`。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_common_records(input_paths, output_path, *, by='id', ignore_case=False, only_positive_strand=False, check_embedded_seqs=False, threads=1, seqkit_path=None)`。
- **输入：** 至少两个文件和新的输出路径；`by` 取 `id`、`name` 或 `sequence`。
- **输出：** `shared_record_count`、`key`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_common_records

for by in ("id", "name", "sequence"):
    result = seqkit_common_records(["x.fa", "other.fa"], str(out / f"{by}.fa"), by=by)
    print(by, result["shared_record_count"])
result = seqkit_common_records(["x.fa", "other.fa"], str(out / "pos.fa"), by="sequence", only_positive_strand=True)
print("positive strand", result["shared_record_count"])
```

- **示例输出：**

```text
id 2
name 1
sequence 2
positive strand 1
```

## 10) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_replace_header()`

- **作用：** 用正则表达式替换序列名称。
- **处理方法：** 对名称运行 `seqkit replace`。替换串可使用 `$1`、`{nr}`（记录序号，`number_width` 补零）、`{fn}` 与 `{kv}`（来自键值文件的值）。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_replace_header(input_path, output_path, pattern, replacement, *, ignore_case=False, key_value_file=None, key_capture_index=1, keep_key=False, keep_untouched=False, missing_key_replacement=None, number_width=1, only_records_matching=None, only_records_matching_regexp=False, threads=1, seqkit_path=None)`。
- **输入：** 序列文件、新的输出路径、`pattern` 与 `replacement`。
- **输出：** `output`、`output_record_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_replace_header

seqkit_replace_header("x.fa", str(out / "renamed.fa"), "^(\\w)", "seq{nr}_$1", number_width=2)
print([l for l in (out / "renamed.fa").read_text().splitlines() if l.startswith(">")])
```

- **示例输出：**

```text
['>seq01_a first', '>seq02_b second', '>seq03_c third']
```

## 11) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_dump_records()`

- **作用：** 把 BAM 文件中每条记录的比对特征导出为表格。
- **处理方法：** 运行 `seqkit bam -T Dump`。需要 `NM` 标签的特征（`Acc`、`Match`、`Mismatch`）要求 BAM 中带有该标签。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_dump_records(bam_path, tsv_path, *, fields=None, seqkit_path=None)`。
- **输入：** BAM 文件、新的 TSV 路径和可选的 `fields` 列表。
- **输出：** `rows`、`row_count`、`fields`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import pysam
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_bam_dump_records

header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
with pysam.AlignmentFile(str(out / "x.bam"), "wb", header=header) as handle:
    for index, (start, nm) in enumerate(((5, 0), (10, 1))):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"read{index}", "ACGT" * 5, 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigarstring = 0, start, 60, "20M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 20)
        read.set_tag("NM", nm)
        handle.write(read)
result = seqkit_bam_dump_records(str(out / "x.bam"), str(out / "dump.tsv"), fields=["Read", "Pos", "Acc", "MapQual"])
for row in result["rows"]:
    print(row)
```

- **示例输出：**

```text
{'Read': 'read0', 'Pos': 5, 'Acc': 100.0, 'MapQual': 60}
{'Read': 'read1', 'Pos': 10, 'Acc': 95.0, 'MapQual': 60}
```

## 12) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_filter_by_alignment_context()`

- **作用：** 按比对起止位置处的参考序列保留 BAM 记录。
- **处理方法：** 运行 `seqkit bam -T AlnContext`：用 Go 正则表达式检验比对两端 `[pos+left_shift, pos+right_shift)` 的参考序列；`stranded` 对反链 reads 取反向互补，`invert` 反转判断。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_filter_by_alignment_context(bam_path, reference_fasta, output_bam, *, regex_start=None, regex_end=None, left_shift=-10, right_shift=10, stranded=False, invert=False, tsv_path=None, seqkit_path=None)`。
- **输入：** BAM、带索引的参考 FASTA、新的输出 BAM，以及 `regex_start` 和/或 `regex_end`。
- **输出：** `output_bam`、`context_rows`（检查过的窗口与匹配标记）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import pysam
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_bam_filter_by_alignment_context

sequence = "ACGTTGCAAGGCTTAACCGGATATCGCGATTTAAACCCGGGTTTACGTACGT"
reference = out / "ref.fa"
reference.write_text(f">chr1\n{sequence}\n")
pysam.faidx(str(reference))
header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": len(sequence)}]}
with pysam.AlignmentFile(str(out / "x.bam"), "wb", header=header) as handle:
    for index, start in enumerate((5, 10)):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"read{index}", sequence[start:start + 20], 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigarstring = 0, start, 60, "20M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 20)
        handle.write(read)
result = seqkit_bam_filter_by_alignment_context(str(out / "x.bam"), str(reference), str(out / "kept.bam"),
                                                regex_start="^GC", left_shift=0, right_shift=3, tsv_path=str(out / "ctx.tsv"))
for row in result["context_rows"]:
    print(row["Read"], row["StartSeq"], row["StartMatch"])
print(pysam.AlignmentFile(str(out / "kept.bam")).count(until_eof=True))
```

- **示例输出：**

```text
read0 GCA -1
read1 GCT -1
2
```

## 13) `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_accuracy_stats()`

- **作用：** BAM 文件的平均比对准确率，含普通平均和按比对长度加权平均。
- **处理方法：** 运行 `seqkit bam -T AccStats`；单条记录的准确率为 `(1 - NM / 比对长度) * 100`。
- **API：** `dnakit.fasta_fastq.seqkit_p15_cli_functions.seqkit_bam_accuracy_stats(bam_path, *, tsv_path, seqkit_path=None)`。
- **输入：** 带 `NM` 标签的 BAM 和新的 TSV 路径。
- **输出：** `accuracy_mean`、`accuracy_weighted_mean`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import pysam
from dnakit.fasta_fastq.seqkit_p15_cli_functions import seqkit_bam_accuracy_stats

header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100}]}
with pysam.AlignmentFile(str(out / "x.bam"), "wb", header=header) as handle:
    for index, (start, nm) in enumerate(((5, 0), (10, 1))):
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = f"read{index}", "ACGT" * 5, 0
        read.reference_id, read.reference_start, read.mapping_quality, read.cigarstring = 0, start, 60, "20M"
        read.query_qualities = pysam.qualitystring_to_array("I" * 20)
        read.set_tag("NM", nm)
        handle.write(read)
result = seqkit_bam_accuracy_stats(str(out / "x.bam"), tsv_path=str(out / "acc.tsv"))
print(result["accuracy_mean"], result["accuracy_weighted_mean"])
```

- **示例输出：**

```text
97.5 97.5
```

## 14) `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_overlap_analysis()`

- **作用：** 通过重叠分析修剪接头、校正错配碱基并保存双端 reads 的重叠部分。
- **处理方法：** 运行 `after.py`，默认开启重叠分析；`analyse_overlap=False` 即 `--no_overlap`，`correct_mismatches=False` 即 `--no_correction`，`mask_mismatch` 即 `--mask_mismatch`，`store_overlap` 即 `--store_overlap on`。
- **API：** `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_overlap_analysis(read1, read2, output_directory, *, correct_mismatches=True, mask_mismatch=False, analyse_overlap=True, store_overlap=False, trim_front=-1, trim_tail=-1, qualified_quality_phred=15, unqualified_base_limit=60, poly_size_limit=35, n_base_limit=5, seq_len_req=35, gzip_output=False, afterqc_path=None)`。
- **输入：** 同一文库的两个 FASTQ 文件和新的输出文件夹。
- **输出：** `good_files`、`bad_files`、`qc_files`、`overlap_files`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.sequencing_qc.afterqc_p15_cli_functions import afterqc_overlap_analysis

result = afterqc_overlap_analysis("sample_R1.fq", "sample_R2.fq", "afterqc_out", store_overlap=True)
print(result["good_files"])
```

## 15) `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_debubble()`

- **作用：** 过滤前去除测序仪气泡区域的 reads。
- **处理方法：** 对 Illumina FASTQ 目录运行 `after.py --debubble`；坏 reads 密度异常的圆形区域内的 reads 被移到 `bad`。
- **API：** `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_debubble(input_directory, output_directory, *, draw=True, trim_front=-1, trim_tail=-1, qualified_quality_phred=15, unqualified_base_limit=60, poly_size_limit=35, n_base_limit=5, seq_len_req=35, gzip_output=False, afterqc_path=None)`。
- **输入：** 输入目录和新的输出文件夹。
- **输出：** `good_files`、`bad_files`、`debubble_files`、`circles`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.sequencing_qc.afterqc_p15_cli_functions import afterqc_debubble

result = afterqc_debubble("lane_fastq_dir", "afterqc_debubble_out", draw=False)
print(result["circles"])
```

## 16) `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_move_barcodes_to_read_names()`

- **作用：** 把单分子条形码从 reads 开头移到序列名中。
- **处理方法：** 运行 `after.py --barcode on`：开头 `barcode_length` 个碱基为条形码，其后的碱基必须等于 `barcode_verify`。文件名含 `barcode_flag` 时 AfterQC 才把它当作条形码文件。
- **API：** `dnakit.sequencing_qc.afterqc_p15_cli_functions.afterqc_move_barcodes_to_read_names(read1, output_directory, *, barcode_length=12, barcode_verify='CAGTA', barcode_flag='barcode', read2=None, qualified_quality_phred=15, unqualified_base_limit=60, poly_size_limit=35, n_base_limit=5, seq_len_req=35, gzip_output=False, afterqc_path=None)`。
- **输入：** 文件名含 `barcode` 的 FASTQ 和新的输出文件夹。
- **输出：** `good_files`、`bad_files`、`qc_files`。
- **示例调用（需要外部程序，文档构建时不执行）：**

```python
from dnakit.sequencing_qc.afterqc_p15_cli_functions import afterqc_move_barcodes_to_read_names

result = afterqc_move_barcodes_to_read_names("run_barcode_R1.fq", "afterqc_barcode_out", barcode_length=12, barcode_verify="CAGTA")
print(result["good_files"])
```

## 17) `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_process_reads()`

- **作用：** 用 Trimmomatic 步骤（滑动窗口、MAXINFO、AVGQUAL、BASECOUNT、MINLEN、MAXLEN、TOPHRED33/64、CROP、HEADCROP、TAILCROP、LEADING、TRAILING）修剪并过滤单端或双端 reads。
- **处理方法：** 按 Trimmomatic 手册建议的顺序组装步骤，以 `-summary` 运行 `SE` 或 `PE`。Trimmomatic 是 Java 程序：传入 `trimmomatic_jar`（`.jar` 文件），或把 `trimmomatic` 启动脚本放到 PATH。
- **API：** `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_process_reads(read1, output_directory, *, read2=None, input_phred=None, convert_quality=None, crop=None, head_crop=None, tail_crop=None, leading=None, trailing=None, sliding_window=None, max_info=None, average_quality=None, base_count=None, min_length=None, max_length=None, threads=1, trimmomatic_jar=None)`。
- **输入：** `read1`（及 `read2`）、新的输出文件夹和至少一个步骤参数；`input_phred` 指定输入质量编码。
- **输出：** `steps`、`survival`（输入、保留和丢弃的 reads 数）、`read_outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import os
from dnakit.sequencing_qc.trimmomatic_p15_cli_functions import trimmomatic_process_reads

jar = os.environ.get("TRIMMOMATIC_JAR")
single = trimmomatic_process_reads("r1.fq", str(out / "se"), sliding_window=[4, 20], min_length=25, max_length=100,
                                   average_quality=20, base_count=["N", 0, 2], input_phred=33, trimmomatic_jar=jar)
print(single["steps"])
print(single["survival"]["input_reads"], single["survival"]["surviving_reads"])
paired = trimmomatic_process_reads("r1.fq", str(out / "pe"), read2="r2.fq", convert_quality=64, input_phred=33, trimmomatic_jar=jar)
print(paired["survival"]["both_surviving_reads"])
```

- **示例输出：**

```text
['SLIDINGWINDOW:4:20', 'AVGQUAL:20', 'BASECOUNT:N:0:2', 'MINLEN:25', 'MAXLEN:100']
2 2
2
```

## 18) `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_pair_reads()`

- **作用：** 按序列名配对 FASTQ 的两个 mate，并拆出没有配对的 reads。
- **处理方法：** 从 jar 中运行 Trimmomatic 的 `Pairomatic` 类。默认比较完整名称，给出 `delimiter` 时只比较其前面的部分（`name/1` 与 `name/2` 用 `/`）；两个文件共有的 reads 必须顺序一致。
- **API：** `dnakit.sequencing_qc.trimmomatic_p15_cli_functions.trimmomatic_pair_reads(read1, read2, output_directory, *, delimiter=None, trimmomatic_jar=None)`。
- **输入：** 两个 FASTQ 文件、新的输出文件夹和 `trimmomatic_jar`（必需）。
- **输出：** 每个输出文件的 `read_counts` 与 `read_outputs`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

out = Path(tempfile.mkdtemp())
import os
from dnakit.sequencing_qc.trimmomatic_p15_cli_functions import trimmomatic_pair_reads

result = trimmomatic_pair_reads("m1.fq", "m2.fq", str(out / "paired"), delimiter="/", trimmomatic_jar=os.environ["TRIMMOMATIC_JAR"])
print(result["read_counts"])
```

- **示例输出：**

```text
{'R1.paired.fastq': 2, 'R1.unpaired.fastq': 1, 'R2.paired.fastq': 2, 'R2.unpaired.fastq': 0}
```
