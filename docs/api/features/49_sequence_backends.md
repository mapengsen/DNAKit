# 可选序列后端：Biopython、scikit-bio、pyfastx、screed 与 py2bit

这一组函数是 `dnakit.sequence` 中基于可选依赖的适配器。EMBL 与 GenBank 记录通过 [Biopython](https://github.com/biopython/biopython) 的 `Bio.SeqIO` 解析、提取和转换；DNA 翻译、k-mer 频率和简并碱基展开使用 [scikit-bio](https://github.com/scikit-bio/scikit-bio)；FASTQ 质量编码通过 [pyfastx](https://github.com/lmdu/pyfastx) 推断；记录数据库通过 [screed](https://github.com/dib-lab/screed) 写入和查询；2bit 文件的区间读取和统计通过 [py2bit](https://github.com/deeptools/py2bit) 完成。所有函数只处理显式给出的文件或序列，不访问网络，也不计算引物、比对或限制性酶切。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

下面的示例直接从 `dnakit.sequence` 导入函数，并使用这两个文件。EMBL 文件只有 45 个碱基和两个 CDS 特征，FASTQ 文件只有两条短读段，它们只用于演示接口。2bit 文件是二进制文件，因此由示例代码自己写出。

**示例文件：** `records.embl`

```text
ID   AB000001; SV 1; linear; genomic DNA; STD; UNC; 45 BP.
XX
DE   test record for dnakit
XX
FH   Key             Location/Qualifiers
FH
FT   CDS             1..18
FT                   /gene="geneA"
FT   CDS             complement(21..32)
FT                   /gene="geneB"
XX
SQ   Sequence 45 BP;
     atggccaaat ttgggcccta aggatccatg aaattttaac cggta                        45
//
```

**示例文件：** `reads.fq`

```text
@q1 desc
ACGTAC
+
IIII#!
@q2
GGCCAA
+
ABCDEF
```

## 1) `dnakit.sequence.read_sequence_records()`

- **作用：** 解析 EMBL 或 GenBank 记录，返回每条记录的序列、长度和特征数量。
- **处理方法：** 用 `Bio.SeqIO.parse` 逐条读取。记录数或碱基数超过 `max_records`、`max_bases` 时直接报错，不截断结果。
- **API：** `dnakit.sequence.read_sequence_records(source, *, format="embl", max_records=1000, max_bases=10000000)`。
- **输入：** EMBL 或 GenBank 文件路径；`format` 取 `"embl"` 或 `"genbank"`。
- **输出：** 含 `record_count`、`total_bases` 和 `records` 的字典。每条记录包括 `id`、`name`、`description`、`length`、`feature_counts` 和 `sequence`。
- **示例代码：**

```python
from dnakit.sequence import read_sequence_records

result = read_sequence_records("records.embl", format="embl")
record = result["records"][0]
print(record["id"], record["length"], record["feature_counts"])
```

- **示例输出：**

```text
AB000001.1 45 {'CDS': 2}
```

## 2) `dnakit.sequence.extract_feature_sequences()`

- **作用：** 提取指定类型特征的序列，例如 CDS 或 gene。
- **处理方法：** 调用 `Bio.SeqFeature` 的 `extract()`。负链特征取反向互补，`join()` 位置按各段拼接。输出的 `start` 和 `end` 覆盖整个特征，为 1-based 闭区间。
- **API：** `dnakit.sequence.extract_feature_sequences(source, *, format="genbank", feature_type="CDS", max_records=1000, max_features=10000, max_bases=10000000)`。
- **输入：** EMBL 或 GenBank 文件路径；`feature_type` 为特征类型名称，不能包含空白。
- **输出：** `feature_count` 和 `features`。每个特征包括 `record_id`、`feature_index`、`gene`、`locus_tag`、`start`、`end`、`strand`、`parts`、`length` 和 `sequence`。
- **示例代码：**

```python
from dnakit.sequence import extract_feature_sequences

result = extract_feature_sequences("records.embl", format="embl", feature_type="CDS")
for feature in result["features"]:
    print(feature["gene"], feature["start"], feature["end"], feature["sequence"])
```

- **示例输出：**

```text
geneA 1 18 ATGGCCAAATTTGGGCCC
geneB 21 32 TTCATGGATCCT
```

## 3) `dnakit.sequence.convert_sequence_file()`

- **作用：** 把 EMBL 或 GenBank 记录转换为 FASTA、EMBL 或 GenBank。
- **处理方法：** 用 `Bio.SeqIO.write` 按记录写出。目标文件必须不存在；超过 `max_records` 时删除未写完的文件。FASTA 输出只保留标题和序列。
- **API：** `dnakit.sequence.convert_sequence_file(source, target, *, input_format="genbank", output_format="fasta", max_records=1000)`。
- **输入：** 源文件路径、新的目标路径、输入格式（`"embl"` 或 `"genbank"`）和输出格式（`"fasta"`、`"embl"` 或 `"genbank"`）。
- **输出：** `record_count`、`source`、`target`、`input_format` 和 `output_format`。
- **示例代码：**

```python
import os

from dnakit.sequence import convert_sequence_file

if os.path.exists("records.fasta"):
    os.remove("records.fasta")  # conversion never overwrites a file
result = convert_sequence_file(
    "records.embl", "records.fasta", input_format="embl", output_format="fasta"
)
with open("records.fasta", encoding="utf-8") as handle:
    print(result["record_count"], handle.readline().strip())
```

- **示例输出：**

```text
1 >AB000001.1 test record for dnakit
```

## 4) `dnakit.sequence.translate_dna_with_genetic_code()`

- **作用：** 用指定的 NCBI 遗传密码表翻译 DNA。
- **处理方法：** 调用 scikit-bio 的 `DNA.translate`，先转录为 RNA 再翻译。`reading_frame` 取 1、2、3（正链）或 -1、-2、-3（反向互补）。`start` 和 `stop` 取 `"ignore"`、`"require"` 或 `"optional"`；`"require"` 找不到对应密码子时报错。含 gap 的序列会被拒绝。
- **API：** `dnakit.sequence.translate_dna_with_genetic_code(sequence, *, genetic_code=1, reading_frame=1, start="ignore", stop="ignore", max_length=1000000)`。
- **输入：** DNA 字符串，允许 IUPAC 简并碱基；`genetic_code` 取 scikit-bio 支持的 NCBI 表编号，例如 1、2 或 11。
- **输出：** `protein`、`protein_length`、`genetic_code`（含 `table_id` 和 `name`）、`reading_frame`、`start` 和 `stop`。
- **示例代码：**

```python
from dnakit.sequence import translate_dna_with_genetic_code

result = translate_dna_with_genetic_code("ATGTGAATAAGA", genetic_code=2)
print(result["protein"], result["genetic_code"]["name"])
```

- **示例输出：**

```text
MWM* Vertebrate Mitochondrial
```

## 5) `dnakit.sequence.calculate_kmer_frequencies()`

- **作用：** 统计 k-mer 的计数或相对频率。
- **处理方法：** 调用 scikit-bio 的 `DNA.kmer_frequencies`。`overlap=False` 时 k-mer 不重叠；`relative=True` 时数值为计数除以参与统计的 k-mer 总数。含 N 等简并碱基的 k-mer 按原样计数。
- **API：** `dnakit.sequence.calculate_kmer_frequencies(sequence, *, k=3, overlap=True, relative=False, max_length=1000000)`。
- **输入：** DNA 字符串；`k` 为 1 到 32 之间的整数。
- **输出：** `kmers`（按序列排序的字典）、`distinct_kmers`、`k`、`overlap`、`relative` 和 `sequence_length`。
- **示例代码：**

```python
from dnakit.sequence import calculate_kmer_frequencies

result = calculate_kmer_frequencies("ACACATTTATTA", k=3, overlap=False, relative=True)
print(result["kmers"])
```

- **示例输出：**

```text
{'ACA': 0.25, 'CAT': 0.25, 'TTA': 0.5}
```

## 6) `dnakit.sequence.expand_degenerate_sequences()`

- **作用：** 列出简并 DNA 序列代表的全部具体序列。
- **处理方法：** 调用 scikit-bio 的 `DNA.expand_degenerates`。展开前先计算结果数，即各简并位点可选碱基数的乘积；超过 `max_sequences` 时直接报错，不生成序列。结果按字典序排列。
- **API：** `dnakit.sequence.expand_degenerate_sequences(sequence, *, max_sequences=1024, max_length=1000000)`。
- **输入：** IUPAC DNA 字符串，例如含 R、Y 或 N 的序列。
- **输出：** `sequences`、`expansion_count`、`degenerate_positions` 和 `sequence_length`。
- **示例代码：**

```python
from dnakit.sequence import expand_degenerate_sequences

print(expand_degenerate_sequences("TRG")["sequences"])
```

- **示例输出：**

```text
['TAG', 'TGG']
```

## 7) `dnakit.sequence.detect_fastq_quality_encoding()`

- **作用：** 推断 FASTQ 的质量偏移（Phred+33 或 Phred+64）以及质量分数范围。
- **处理方法：** 调用 pyfastx 的 `Fastq`，不建立索引，扫描整个文件后取其报告的偏移和字符范围。pyfastx 返回 0 时表示 33 和 64 都与观察到的字符兼容，此时 `phred_offset` 为 `None`，且不给出分数范围。扫描前先检查文件大小是否超过 `max_bytes`。
- **API：** `dnakit.sequence.detect_fastq_quality_encoding(source, *, max_bytes=1000000000)`。
- **输入：** FASTQ 文件路径，支持 gzip 压缩文件。
- **输出：** `phred_offset`、`possible_phred_offsets`、`min_quality_ascii`、`max_quality_ascii`、`min_quality_score`、`max_quality_score`、`min_length` 和 `max_length`。
- **示例代码：**

```python
from dnakit.sequence import detect_fastq_quality_encoding

result = detect_fastq_quality_encoding("reads.fq")
print(result["phred_offset"], result["min_quality_score"], result["max_quality_score"])
```

- **示例输出：**

```text
33 0 40
```

## 8) `dnakit.sequence.build_screed_database()`

- **作用：** 把 FASTA 或 FASTQ 文件写成新的 screed SQLite 记录数据库。
- **处理方法：** 调用 `screed.create_db`，字段定义取自 screed 的 FASTA 或 FASTQ 读取器。若路径不以 `_screed` 结尾，会自动补上该后缀。目标已存在时报错；超过 `max_records` 时删除未完成的数据库。
- **API：** `dnakit.sequence.build_screed_database(source, *, database, max_records=1000)`。
- **输入：** FASTA 或 FASTQ 文件路径，以及一个尚不存在的数据库路径。
- **输出：** `database`（实际写出的路径）、`format`（`"fasta"` 或 `"fastq"`）和 `record_count`。
- **示例代码：**

```python
import os

from dnakit.sequence import build_screed_database

if os.path.exists("reads_screed"):
    os.remove("reads_screed")  # building never overwrites a database
result = build_screed_database("reads.fq", database="reads")
print(result["format"], result["record_count"])
```

- **示例输出：**

```text
fastq 2
```

## 9) `dnakit.sequence.fetch_screed_records()`

- **作用：** 按名称读取 screed 数据库中的记录；不指定名称时读取全部记录。
- **处理方法：** 用 `screed.ScreedDB` 按键查询。序列统一转为大写。FASTQ 记录额外返回 `quality`。缺失的名称会报错，并列出前 10 个缺失名称。
- **API：** `dnakit.sequence.fetch_screed_records(database, *, names=None, max_records=1000)`。
- **输入：** `build_screed_database()` 返回的数据库路径，或不带 `_screed` 后缀的同一路径；`names` 为记录名列表。
- **输出：** `record_count` 和 `records`。每条记录包括 `name`、`description`、`length` 和 `sequence`；FASTQ 记录还包括 `quality`。
- **示例代码：**

```python
from dnakit.sequence import fetch_screed_records

record = fetch_screed_records("reads", names=["q1"])["records"][0]
print(record["name"], record["sequence"], record["quality"])
```

- **示例输出：**

```text
q1 ACGTAC IIII#!
```

## 10) `dnakit.sequence.summarize_twobit()`

- **作用：** 统计 2bit 文件中每条序列的长度、硬屏蔽（N）碱基数、ACGT 碱基数和 GC 含量。
- **处理方法：** 调用 py2bit 的 `chroms()`、`bases()` 和 `hardMaskedBlocks()`。软屏蔽（小写）被忽略。GC 百分比以 ACGT 碱基为分母，N 不计入。
- **API：** `dnakit.sequence.summarize_twobit(source, *, max_chromosomes=10000)`。
- **输入：** 2bit 文件路径。下面的示例先写出一个小文件 `toy.2bit`，只包含 ACGT 碱基。
- **输出：** `chromosome_count`、`total_length`、`total_hard_masked` 和 `chromosomes`。每条序列包括 `name`、`length`、`acgt_bases`、`hard_masked_bases` 和 `gc_percent_of_acgt`。
- **示例代码：**

```python
import struct

from dnakit.sequence import summarize_twobit

# Minimal little-endian .2bit writer for this example: ACGT only, no masks.
records = [("chr1", "ACGTACGTGGCC"), ("chr2", "TTTTAAAACCCCGGGG")]
codes = {"T": 0, "C": 1, "A": 2, "G": 3}
index, body = b"", b""
offset = 16 + sum(1 + len(name) + 4 for name, _ in records)
for name, sequence in records:
    packed = bytearray((len(sequence) + 3) // 4)
    for position, base in enumerate(sequence):
        packed[position // 4] |= codes[base] << (6 - 2 * (position % 4))
    block = struct.pack("<IIII", len(sequence), 0, 0, 0) + bytes(packed)
    index += bytes([len(name)]) + name.encode() + struct.pack("<I", offset)
    body += block
    offset += len(block)
with open("toy.2bit", "wb") as handle:
    handle.write(struct.pack("<IIII", 0x1A412743, 0, len(records), 0) + index + body)

result = summarize_twobit("toy.2bit")
print([(c["name"], c["length"], c["gc_percent_of_acgt"]) for c in result["chromosomes"]])
```

- **示例输出：**

```text
[('chr1', 12, 66.66666666666667), ('chr2', 16, 50.0)]
```

## 11) `dnakit.sequence.fetch_twobit_region()`

- **作用：** 读取 2bit 文件中某条序列的一个区间。
- **处理方法：** 调用 py2bit 的 `sequence()`。坐标为 0-based、左闭右开。`end` 超过序列长度时报错，而不是像 py2bit 那样静默截断。输出统一为大写。
- **API：** `dnakit.sequence.fetch_twobit_region(source, chromosome, start, end, *, max_length=1000000)`。
- **输入：** 2bit 文件路径、序列名、0-based 的 `start` 和左闭右开的 `end`。示例读取上一节写出的 `toy.2bit`。
- **输出：** `chromosome`、`start`、`end`、`length`、`hard_masked_bases` 和 `sequence`。
- **示例代码：**

```python
from dnakit.sequence import fetch_twobit_region

print(fetch_twobit_region("toy.2bit", "chr1", 2, 8)["sequence"])
```

- **示例输出：**

```text
GTACGT
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-formats]"     # biopython：EMBL/GenBank 的读取、提取和转换
python -m pip install "dnakit[bio-skibio]"      # scikit-bio：翻译、k-mer 与简并碱基展开
python -m pip install "dnakit[bio-pyfastx]"     # pyfastx：FASTQ 质量编码
python -m pip install "dnakit[bio-screed]"      # screed：SQLite 记录数据库
python -m pip install "dnakit[bio-py2bit]"      # py2bit：2bit 统计与区间读取
```

各包的项目地址：[Biopython](https://github.com/biopython/biopython)、[scikit-bio](https://github.com/scikit-bio/scikit-bio)、[pyfastx](https://github.com/lmdu/pyfastx)、[screed](https://github.com/dib-lab/screed)、[py2bit](https://github.com/deeptools/py2bit)。

## 验证方式

测试文件是 `others/tests/unit/test_sequence_backends.py`，共 32 项。主要检查：

- EMBL/GenBank 的读取、特征提取和格式转换：基因序列、负链反向互补和 `join()` 位置都与手工计算的结果一致。输入文件由 Biopython 在测试中生成。
- scikit-bio 包装：不同遗传密码表和读框的翻译、起始与终止策略、非重叠和相对频率的 k-mer 计数，以及简并碱基展开的数量和内容，都与手工计算一致。重叠 k-mer 计数还与直接调用 scikit-bio 的结果逐项相同。
- pyfastx、screed 和 py2bit：FASTQ 质量偏移与分数范围、screed 的 FASTA/FASTQ 写入与读取、2bit 的 GC 与硬屏蔽计数，以及 0-based 半开区间读取，都与手工计算一致。
- 超过 `max_records`、`max_bases`、`max_features`、`max_length`、`max_sequences` 或 `max_bytes` 时报错；转换或建库失败时不留下未完成的输出文件。
- 缺少可选依赖时给出安装提示；gap 字符、未知遗传密码表编号、非 2bit 文件、非 FASTQ 输入和已存在的输出路径都会被拒绝。

## 当前范围

这一批只包装在 dnakit 中尚未以相同方式提供的 DNA 和基因特征操作。以下操作没有包装：

- pyfaidx 的区间读取已由 `dnakit.integrations.fetch_pyfaidx_region()` 提供。
- pyfastx 的 FASTA/FASTQ 汇总统计已由 `dnakit.integrations.summarize_pyfastx()` 提供；dnaio 的流式汇总已由 `dnakit.integrations.summarize_dnaio_fastx()` 提供。dnaio 的写入、配对和反向互补由原生的 `dnakit.bioinformatics.reads` 与 `dnakit.fasta_fastq` 提供。
- FASTA/FASTQ 解析（Biopython `Bio.SeqIO`）与原生 `dnakit.bioinformatics.reads.read_fastx()` 重叠，未包装。GenBank 解析与原生 `dnakit.io.genbank` 重叠。这里保留它，是因为特征提取和格式转换与 Biopython 的 EMBL 解析共用同一条路径。
- 反向互补和 GC 含量由原生 `complement_dna`、`reverse_complement_dna` 和 `calculate_sequence_metrics` 提供，未包装 scikit-bio 的同类方法。
- 原生 `translate_dna_to_protein()` 只支持标准遗传密码表，本批的翻译函数用 scikit-bio 支持多种密码表。Biopython 的 `Seq.translate` 与本批的翻译功能重复，未包装。
- 熔解温度（Biopython `MeltingTemp`）与原生 `dnakit.thermodynamics.melting_temperature()` 功能重叠，未包装。
- 汉明距离由原生 `dnakit.similarity.hamming_distance()` 提供，未包装 scikit-bio 的 `hamming`。
- 限制性酶切由原生 `dnakit.patterns.scan_restriction_sites()` 和 `dnakit.molbio.digest_restriction()` 提供，未包装 Biopython 的 `Bio.Restriction`。
- GC 偏移（GC skew）已由原生 descriptors 提供，未包装。
- py2bit 只能读取，没有写入接口，因此本批没有 2bit 写入函数。
- pyfastx 的随机访问需要在输入文件旁边建立 `.fxi` 索引，本批不使用这一功能。
