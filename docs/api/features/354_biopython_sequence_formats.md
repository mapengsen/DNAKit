# Biopython 序列格式与 Seq 辅助函数：QUAL、IMGT、GFA、nib、PIR、SeqXML、SFF、Xdna、Entrez XML 与 BGZF 偏移

这一组函数位于 `dnakit.sequence.biopython_p18_functions`，封装 [Biopython](https://biopython.org/)（`Bio.SeqIO`、`Bio.Seq`、`Bio.SeqFeature`、`Bio.Entrez`、`Bio.bgzf`、`Bio.UniProt.GOA`、`Bio.Emboss`、`Bio.Align.substitution_matrices`），补齐此前 Biopython 封装未覆盖的操作：读写较少见的序列文件格式、对序列字符串的精确匹配操作、INSDC 位置字符串、EMBOSS `primersearch` 输入、GO 注释过滤、替换矩阵文本、已保存的 NCBI Entrez XML 文件以及 BGZF 虚拟偏移。

读取函数只读本地文件。写入函数从不覆盖：输出路径必须是新文件。序列字符串函数保持输入的大小写与字符，与 `Bio.Seq` 一致。所有函数都不联网；`parse_entrez_xml()` 遇到本地缺失的 DTD 会报错而不是下载。Biopython 是可选依赖（`pip install "dnakit[bio-formats]"`），只有调用函数时才导入；缺少时报错并给出安装命令。示例从模块 `dnakit.sequence.biopython_p18_functions` 导入。

## 1) `dnakit.sequence.biopython_p18_functions.read_sequence_file()`

- **作用：** 读取较少见的 `Bio.SeqIO` 格式，并把记录作为普通字典返回。
- **处理方法：** 按所选格式调用 `Bio.SeqIO.parse()`：`qual`（只有 PHRED 质量，序列未定义，返回 null）、`imgt`、`genbank-cds` 与 `embl-cds`（每个 CDS 特征一条记录，内容为翻译产物）、`gck`、`gfa1`、`gfa2`（组装图片段）、`ig`、`nib`、`pir`、`seqxml`、`tab`、`xdna`、`abi` 与 `abi-trim`（Biopython 的 Mott 式修剪，错误概率阈值 0.05）。ABI 文件的原始 ABIF 数据块默认不返回，需要 `include_raw_trace=True`。文件记录数超过 `max_records` 时报错。
- **API：** `dnakit.sequence.biopython_p18_functions.read_sequence_file(path, *, file_format='tab', max_records=10000, include_sequence=True, include_raw_trace=False)`。
- **输入：** 本地文件路径、`file_format`，以及可选的 `max_records`（默认 10000）、`include_sequence` 和 `include_raw_trace`。
- **输出：** `record_count` 与 `records`；每条记录含 `id`、`name`、`description`、`length`、`sequence`、`annotations`、`letter_annotations`（QUAL 为 `phred_quality` 列表）、`feature_count`、`features`（类型、位置、限定符）和 `dbxrefs`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import read_sequence_file

Path("records.tab").write_text("seq1\tACGTACGTAC\nseq2\tGGTTAACC\n", encoding="utf-8")
result = read_sequence_file("records.tab", file_format="tab")
print(result["record_count"], [(r["id"], r["length"]) for r in result["records"]])
print(result["records"][0]["sequence"])
```

- **示例输出：**

```text
2 [('seq1', 10), ('seq2', 8)]
ACGTACGTAC
```

## 2) `dnakit.sequence.biopython_p18_functions.write_sequence_file()`

- **作用：** 把序列记录写成新的 IMGT、nib、PHD、PIR、SeqXML、tab 或 Xdna 文件。
- **处理方法：** 由映射构造 `SeqRecord` 对象并调用 `Bio.SeqIO.write()`。`nib` 与 `xdna` 格式只能保存一条记录。SFF 文件含流信号数据，因此由 `convert_sequence_format()` 从 SFF 源文件写出。输出路径必须不存在；写入失败会删除不完整文件。
- **API：** `dnakit.sequence.biopython_p18_functions.write_sequence_file(records, output_path, *, file_format='tab')`。
- **输入：** 记录映射列表（含 `id` 与 `sequence`，可选 `name`、`description`、`molecule_type`（`DNA` 或 `RNA`）、`topology`（`linear` 或 `circular`）、`qualities`（PHRED 整数，`phd` 必需））、新的输出路径和 `file_format`。
- **输出：** `output_path`、`file_format`、`record_count` 和 `bytes_written`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import write_sequence_file

out = Path(tempfile.mkdtemp()) / "demo.pir"
records = [
    {"id": "a", "sequence": "ACGTACGT", "description": "first"},
    {"id": "b", "sequence": "GGTTAACC"},
]
result = write_sequence_file(records, out, file_format="pir")
print(result["record_count"], result["bytes_written"])
print(out.read_text(encoding="utf-8"))
```

- **示例输出：**

```text
2 44
>D1;a
a - first
ACGTACGT*
>D1;b
b
GGTTAACC*
```

## 3) `dnakit.sequence.biopython_p18_functions.convert_sequence_format()`

- **作用：** 在 `Bio.SeqIO` 格式之间转换序列文件，包括输出 SFF。
- **处理方法：** 用 `Bio.SeqIO.parse()` 解析源文件，用 `Bio.SeqIO.write()` 写出。源格式为 `read_sequence_file()` 的格式再加 embl、fasta、fastq、genbank、sff 和 sff-trim；目标格式为 embl、fasta、fastq、genbank、imgt、nib、phd、pir、qual、seqxml、sff、tab 和 xdna。SFF 目标使用 `SffWriter`，并复制 SFF 源文件中的 Roche XML 清单。目标格式无法保存的信息（FASTA 中的特征、tab 中的质量）会丢失，与 Biopython 一致。
- **API：** `dnakit.sequence.biopython_p18_functions.convert_sequence_format(source, output_path, *, source_format='genbank', target_format='tab', max_records=10000)`。
- **输入：** 源路径、新的输出路径、`source_format`、`target_format` 与 `max_records`。
- **输出：** `record_count`、`source_format`、`target_format`、`roche_manifest_copied` 以及路径。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import convert_sequence_format

Path("records.tab").write_text("seq1\tACGTACGTAC\nseq2\tGGTTAACC\n", encoding="utf-8")
out = Path(tempfile.mkdtemp()) / "records.fasta"
result = convert_sequence_format("records.tab", out, source_format="tab", target_format="fasta")
print(result["record_count"], result["target_format"])
print(out.read_text(encoding="utf-8"))
```

- **示例输出：**

```text
2 fasta
>seq1
ACGTACGTAC
>seq2
GGTTAACC
```

## 4) `dnakit.sequence.biopython_p18_functions.read_sff_manifest()`

- **作用：** 返回 SFF 文件中保存的 Roche XML 清单。
- **处理方法：** 对二进制文件调用 `Bio.SeqIO.SffIO.ReadRocheXmlManifest()`。没有清单的文件返回 `has_manifest = False`；不是 SFF 的文件报错。
- **API：** `dnakit.sequence.biopython_p18_functions.read_sff_manifest(path)`。
- **输入：** 本地 SFF 文件路径。
- **输出：** `has_manifest`、`xml` 和 `xml_length`。
- **示例代码：**

```python
from Bio.Seq import Seq
from Bio.SeqIO import SffIO
from Bio.SeqRecord import SeqRecord

from dnakit.sequence.biopython_p18_functions import read_sff_manifest

record = SeqRecord(Seq("ACGT"), id="read1")
record.letter_annotations["phred_quality"] = [30, 31, 32, 33]
record.annotations.update(
    flow_values=(0, 100, 100, 100, 100, 0, 0, 0), flow_index=(2, 1, 1, 1), flow_chars="TACGTACG",
    flow_key="", clip_qual_left=1, clip_qual_right=4, clip_adapter_left=0, clip_adapter_right=0,
)
with open("demo.sff", "wb") as handle:
    SffIO.SffWriter(handle, xml="<manifest><run/></manifest>").write_file([record])
result = read_sff_manifest("demo.sff")
print(result["has_manifest"], result["xml"], result["xml_length"])
```

- **示例输出：**

```text
True <manifest><run/></manifest> 27
```

## 5) `dnakit.sequence.biopython_p18_functions.count_subsequence()`

- **作用：** 统计精确子序列的出现次数，可选不重叠或重叠计数。
- **处理方法：** `overlap=False` 调用 `Seq.count()`（从左到右不重叠扫描）；`overlap=True` 调用 `Seq.count_overlap()`。匹配区分大小写。`start` 与 `end` 是可选的 0 起始切片边界。
- **API：** `dnakit.sequence.biopython_p18_functions.count_subsequence(sequence, subsequence, *, overlap=False, start=None, end=None)`。
- **输入：** `sequence`、`subsequence`、`overlap`、`start`、`end`。
- **输出：** `count`、`overlap`、`sequence_length` 和 `subsequence`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import count_subsequence

sequence = "AAAAGAAAAGAAAA"
print(count_subsequence(sequence, "AA")["count"])
print(count_subsequence(sequence, "AA", overlap=True)["count"])
print(count_subsequence(sequence, "AA", overlap=True, start=2, end=8)["count"])
```

- **示例输出：**

```text
6
9
3
```

## 6) `dnakit.sequence.biopython_p18_functions.search_subsequences()`

- **作用：** 找出多个精确子序列的所有位置。
- **处理方法：** 调用 `Seq.search()`，报告每个子序列的每个位置；匹配之间可以重叠。位置从 0 开始并按升序排列；匹配区分大小写。
- **API：** `dnakit.sequence.biopython_p18_functions.search_subsequences(sequence, subsequences)`。
- **输入：** `sequence` 以及最多 1000 条 `subsequences`。
- **输出：** `hit_count` 与 `hits`，每个命中含 `position` 和 `subsequence`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import search_subsequences

result = search_subsequences("AAAAGAAAAGAAAA", ["AAG", "GA"])
print(result["hit_count"])
for hit in result["hits"]:
    print(hit["position"], hit["subsequence"])
```

- **示例输出：**

```text
4
2 AAG
4 GA
7 AAG
9 GA
```

## 7) `dnakit.sequence.biopython_p18_functions.split_sequence()`

- **作用：** 按分隔符切分序列。
- **处理方法：** 调用 `Seq.split()`；`from_right=True` 时调用 `Seq.rsplit()`。不给分隔符时按空白切分。`maxsplit` 限制切分次数（-1 表示不限）。
- **API：** `dnakit.sequence.biopython_p18_functions.split_sequence(sequence, separator=None, *, maxsplit=-1, from_right=False)`。
- **输入：** `sequence`，可选的 `separator`、`maxsplit`、`from_right`。
- **输出：** `piece_count` 与 `pieces`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import split_sequence

sequence = "AAAAGAAAAGAAAA"
print(split_sequence(sequence, "G")["pieces"])
print(split_sequence(sequence, "G", maxsplit=1)["pieces"])
print(split_sequence(sequence, "G", maxsplit=1, from_right=True)["pieces"])
```

- **示例输出：**

```text
['AAAA', 'AAAA', 'AAAA']
['AAAA', 'AAAAGAAAA']
['AAAAGAAAA', 'AAAA']
```

## 8) `dnakit.sequence.biopython_p18_functions.trim_sequence()`

- **作用：** 按 `Bio.Seq.Seq` 的方式修剪序列两端。
- **处理方法：** `mode` 为 `strip`（两端）、`lstrip` 或 `rstrip` 时移除 `pattern` 中的任意字符（为 null 时移除空白）；为 `removeprefix` / `removesuffix` 时移除一次完全匹配的 `pattern` 字符串。保持大小写。
- **API：** `dnakit.sequence.biopython_p18_functions.trim_sequence(sequence, *, mode='strip', pattern=None)`。
- **输入：** `sequence`、`mode` 和 `pattern`。
- **输出：** `trimmed`、`trimmed_length`、`original_length`、`mode` 和 `pattern`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import trim_sequence

for mode, pattern in [("strip", "A"), ("lstrip", "A"), ("removeprefix", "AAAA"), ("removesuffix", "AAAA")]:
    print(mode, trim_sequence("AAAAGAAAAGAAAA", mode=mode, pattern=pattern)["trimmed"])
```

- **示例输出：**

```text
strip GAAAAG
lstrip GAAAAGAAAA
removeprefix GAAAAGAAAA
removesuffix AAAAGAAAAG
```

## 9) `dnakit.sequence.biopython_p18_functions.parse_feature_location()`

- **作用：** 解析 INSDC 特征位置字符串，例如 `complement(join(1..10,20..30))`。
- **处理方法：** 调用 `Bio.SeqFeature.Location.fromstring()`。整体起点为 0 起始、终点为开区间；输入字符串则是 1 起始。模糊位置（`<1`、`>200`、`(5.9)`、`12^13`）按 Biopython 的位置类型报告。
- **API：** `dnakit.sequence.biopython_p18_functions.parse_feature_location(location)`。
- **输入：** 位置字符串。
- **输出：** `location_type`、`start` 与 `end`（类型、文本、数值）、`strand`、`operator`（`join`、`order` 或 `single`）、`ref`、`part_count`、`parts` 和规范化文本。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import parse_feature_location

joined = parse_feature_location("complement(join(1..10,20..30))")
print(joined["start"]["value"], joined["end"]["value"], joined["strand"], joined["operator"])
print([(p["start"]["value"], p["end"]["value"]) for p in joined["parts"]])
fuzzy = parse_feature_location("<1..>200")
print(fuzzy["start"]["type"], fuzzy["end"]["type"], fuzzy["normalized"])
print(parse_feature_location("AB000123.1:5..10")["ref"])
```

- **示例输出：**

```text
0 30 -1 join
[(19, 30), (0, 10)]
BeforePosition AfterPosition [<0:>200](+)
AB000123.1
```

## 10) `dnakit.sequence.biopython_p18_functions.format_primersearch_input()`

- **作用：** 由命名的正反向引物生成 EMBOSS `primersearch` 的输入文件。
- **处理方法：** 用 `Bio.Emboss.PrimerSearch.InputRecord.add_primer_set()` 逐对添加，返回文本，每组一行：名称、正向引物、反向引物。引物必须是 IUPAC DNA。给出 `output_path` 时，文本同时写入该新文件。
- **API：** `dnakit.sequence.biopython_p18_functions.format_primersearch_input(primer_sets, *, output_path=None)`。
- **输入：** 引物组名到 `[forward, reverse]` 的映射，以及可选的新 `output_path`。
- **输出：** `primer_set_count`、`text` 和 `output_path`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import format_primersearch_input

result = format_primersearch_input({"p1": ["ACGTAC", "GGTACC"], "p2": ["AAAACG", "TTTTGC"]})
print(result["primer_set_count"])
print(result["text"])
```

- **示例输出：**

```text
2
p1 ACGTAC GGTACC
p2 AAAACG TTTTGC
```

## 11) `dnakit.sequence.biopython_p18_functions.filter_go_annotations()`

- **作用：** 按字段取值保留 GAF 注释文件中的记录。
- **处理方法：** 用 `Bio.UniProt.GOA.gafiterator()` 读取文件，只要所列字段中任意一个取到任一可接受的值就保留该记录（`Bio.UniProt.GOA.record_has()`）。`DB:Reference`、`Taxon_ID` 等列表型字段在任一元素匹配时即匹配。字段名沿用 Biopython 的 GAF 记录：`DB`、`DB_Object_ID`、`DB_Object_Symbol`、`Qualifier`、`GO_ID`、`DB:Reference`、`Evidence`、`With`、`Aspect`、`DB_Object_Name`、`Synonym`、`DB_Object_Type`、`Taxon_ID`、`Date`、`Assigned_By` 等。
- **API：** `dnakit.sequence.biopython_p18_functions.filter_go_annotations(path, field_values, *, max_records=10000)`。
- **输入：** GAF 文件路径、字段名到可接受取值列表的映射，以及 `max_records`。
- **输出：** `records_read`、`kept_count` 与保留的 `records`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import filter_go_annotations

row1 = ["UniProtKB", "P12345", "GENE1", "", "GO:0008150", "PMID:1", "IDA", "", "P", "gene one",
        "", "protein", "taxon:9606", "20200101", "UniProt", "", ""]
row2 = ["UniProtKB", "Q99999", "GENE2", "NOT", "GO:0003674", "PMID:2|PMID:3", "IEA", "", "F", "gene two",
        "", "protein", "taxon:10090", "20200102", "UniProt", "", ""]
Path("demo.gaf").write_text("!gaf-version: 2.1\n" + "\n".join("\t".join(r) for r in (row1, row2)) + "\n", encoding="utf-8")
result = filter_go_annotations("demo.gaf", {"DB:Reference": ["PMID:3"]})
print(result["records_read"], result["kept_count"], result["records"][0]["DB_Object_Symbol"])
print(filter_go_annotations("demo.gaf", {"Aspect": ["P"]})["records"][0]["GO_ID"])
```

- **示例输出：**

```text
2 1 GENE2
GO:0008150
```

## 12) `dnakit.sequence.biopython_p18_functions.format_substitution_matrix()`

- **作用：** 把内置替换矩阵（例如 `NUC.4.4`）返回为格式化文本。
- **处理方法：** 用 `Bio.Align.substitution_matrices.load()` 载入矩阵，再用 `Array.format()` 格式化。`fmt` 是 printf 风格的数字格式（`%i`、`%.1f` 等）。文本先是矩阵头部注释，然后是字母表和每个字母一行。
- **API：** `dnakit.sequence.biopython_p18_functions.format_substitution_matrix(name='NUC.4.4', *, fmt='%i')`。
- **输入：** 内置矩阵 `name` 和 `fmt`。
- **输出：** `name`、`alphabet`、`fmt` 和 `text`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import format_substitution_matrix

result = format_substitution_matrix("NUC.4.4", fmt="%i")
print(result["alphabet"])
print("\n".join(result["text"].splitlines()[-16:-12]))
```

- **示例输出：**

```text
ATGCSWRYKMBVHDN
   A  T  G  C  S  W  R  Y  K  M  B  V  H  D  N
A  5 -4 -4 -4 -4  1  1 -4 -4  1 -4 -1 -1 -1 -2
T -4  5 -4 -4 -4  1 -4  1  1 -4 -1 -4 -1 -1 -2
G -4 -4  5 -4  1 -4  1 -4  1 -4 -1 -1 -4 -1 -2
```

## 13) `dnakit.sequence.biopython_p18_functions.parse_entrez_xml()`

- **作用：** 把已保存的 NCBI Entrez（E-utilities）XML 结果文件解析为普通 Python 数据。
- **处理方法：** `multiple=False` 调用 `Bio.Entrez.read()`（单个结果，如 ESearch、EInfo、ELink）；`multiple=True` 调用 `Bio.Entrez.parse()`（记录列表，如 EFetch 或 ESummary 集合）。Biopython 需要文件中声明的 DTD。本函数不会下载：先读取 DOCTYPE，若 DTD 既不在 Biopython 自带目录也不在用户 DTD 目录，就报错停止。
- **API：** `dnakit.sequence.biopython_p18_functions.parse_entrez_xml(path, *, multiple=False, max_records=1000)`。
- **输入：** 本地 XML 文件路径、`multiple` 与 `max_records`。
- **输出：** `dtd`、`multiple`、`record_count` 和 `data`（嵌套列表与字典）。
- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence.biopython_p18_functions import parse_entrez_xml

Path("esearch.xml").write_text(
    '<?xml version="1.0" encoding="UTF-8" ?>\n'
    '<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" '
    '"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">\n'
    "<eSearchResult><Count>2</Count><RetMax>2</RetMax><RetStart>0</RetStart>"
    "<IdList><Id>111</Id><Id>222</Id></IdList><TranslationSet></TranslationSet>"
    "<QueryTranslation>BRCA1[Gene]</QueryTranslation></eSearchResult>\n",
    encoding="utf-8",
)
result = parse_entrez_xml("esearch.xml")
print(result["dtd"], result["data"]["Count"], result["data"]["IdList"], result["data"]["QueryTranslation"])
```

- **示例输出：**

```text
esearch.dtd 2 ['111', '222'] BRCA1[Gene]
```

## 14) `dnakit.sequence.biopython_p18_functions.bgzf_virtual_offset()`

- **作用：** 合成或拆分 BGZF 虚拟文件偏移。
- **处理方法：** 给出 `block_start`（BGZF 块的压缩偏移，小于 2**48）和 `within_block`（解压块内偏移，小于 2**16），计算 64 位虚拟偏移（`Bio.bgzf.make_virtual_offset()`）；或给出 `virtual_offset`，拆分为两部分（`Bio.bgzf.split_virtual_offset()`）。两种形式必须且只能用一种。
- **API：** `dnakit.sequence.biopython_p18_functions.bgzf_virtual_offset(block_start=None, within_block=None, *, virtual_offset=None)`。
- **输入：** `block_start` 与 `within_block`，或 `virtual_offset`。
- **输出：** `direction`（`make` 或 `split`）、`block_start`、`within_block` 和 `virtual_offset`。
- **示例代码：**

```python
from dnakit.sequence.biopython_p18_functions import bgzf_virtual_offset

made = bgzf_virtual_offset(100, 50)
print(made["direction"], made["virtual_offset"])
split = bgzf_virtual_offset(virtual_offset=made["virtual_offset"])
print(split["direction"], split["block_start"], split["within_block"])
```

- **示例输出：**

```text
make 6553650
split 100 50
```
