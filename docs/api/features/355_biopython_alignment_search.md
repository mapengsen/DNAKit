# Biopython 比对格式、MAF 索引、HMMER 与 Infernal 报告及比对分析

这一组函数位于 `dnakit.alignment.biopython_p18_alignment_functions`，封装 [Biopython](https://biopython.org/)（`Bio.Align`、`Bio.AlignIO.MafIO`、`Bio.SearchIO`、`Bio.Nexus`），补齐此前 Biopython 封装未覆盖的比对操作。它们读写 A2M、bigBed、bigMaf、bigPsl、chain、Clustal、EMBOSS、Exonerate、MAF、Mauve、MSF、NEXUS、PHYLIP、PSL、Stockholm 和 tabular 比对格式；按参考区间查询 MAF 文件；读取 HMMER 与 Infernal 报告，并对检索结果建立索引、过滤和排序；运行 FOGSAA 与密码子感知的成对比对器；分析比对（行排序、替换计数、列频率、一致性计数、坐标映射、合并、映射）以及 DNA 矩阵（恒定列与全空位列、位点状态、子集、插入空位列）。

比对以“序列名到带空位序列（`-` 表示空位）”的映射传入并以同样方式返回，因此所有结果都可序列化为 JSON。写入函数从不覆盖：输出路径必须是新文件。Biopython 是可选依赖（`pip install "dnakit[bio-formats]"`），只有调用函数时才导入。示例从模块 `dnakit.alignment.biopython_p18_alignment_functions` 导入。

有意不包含：`QueryResult.hit_map`、`hsp_map` 与 `Hit.map` 会对每个命中或 HSP 应用任意 Python 函数，而函数对象无法通过 JSON 类型的接口传入，因此该操作被跳过；声明式的情形由 `filter_search_results()` 和 `sort_search_results()` 覆盖。

## 示例数据

示例声明了三个输入：一个在参考序列 `hg.chr1` 上含两个比对块的 MAF 文件、一个只含两个查询的 HMMER `--domtblout` 简化表，以及一个含三个查询的 Infernal `cmscan --tblout` 表（Biopython 无法为 HMMER 表格格式建立索引，所以索引示例使用 Infernal 文件）。

**示例文件：** `demo.maf`

```text
##maf version=1 scoring=none
a score=100
s hg.chr1 100 10 + 1000 ACGTACGTAC
s mm.chr2 200  9 + 2000 ACGTA-GTAC

a score=80
s hg.chr1 120  8 + 1000 GGTTAACC
s mm.chr2 400  8 + 2000 GGTTAACC
```

**示例文件：** `hits.domtab`

```text
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
Globin               PF00042.17   110 gi|4885477|ref|NP_005359.1| -            154     1e-22   80.5   0.3   1   1   1.1e-26   1.6e-22   79.8   0.3     1   109     7   112     7   113 0.97 Globin
Xpo1                 PF08389.7    148 gi|22748937|ref|NP_065801.1| -           1204   8.5e-34  116.6   7.8   1   2   2.3e-37   1.2e-33  116.1   4.9     2   148   110   271   109   271 0.98 Exportin 1-like protein
Xpo1                 PF08389.7    148 gi|22748937|ref|NP_065801.1| -           1204   8.5e-34  116.6   7.8   2   2      0.52   2.6e+03   -1.8   0.0   112   139   499   525   496   529 0.86 Exportin 1-like protein
IBN_N                PF03810.14    77 gi|22748937|ref|NP_065801.1| -           1204    0.0044   16.9   0.0   1   2   7.3e-06     0.036   13.9   0.0     4    75    36    98    33   100 0.87 Importin-beta N-terminal domain
IBN_N                PF03810.14    77 gi|22748937|ref|NP_065801.1| -           1204    0.0044   16.9   0.0   2   2         2     1e+04   -3.5   0.0    57    74   168   185   165   186 0.84 Importin-beta N-terminal domain
Rac1                 PF09632.5    297 gi|22748937|ref|NP_065801.1| -           1204     0.095   11.7   0.3   1   2   0.00026       1.3    8.0   0.3    94   171    89   169    82   254 0.75 Rac1-binding domain
Rac1                 PF09632.5    297 gi|22748937|ref|NP_065801.1| -           1204     0.095   11.7   0.3   2   2     0.041     2e+02    0.8   0.0   150   181   310   341   296   354 0.82 Rac1-binding domain
#
# [ok]
```

**示例文件：** `cm.tbl`

```text
#target name         accession query name              accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#------------------- --------- ----------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
U2                   RF00004   ENA|BK006935|BK006935.2 -          cm        1      193    52929    53083      +    no    1 0.44   0.0   13.5      0.91 ?   U2 spliceosomal RNA
U2                   RF00004   ENA|BK006935|BK006935.2 -          cm        1      193   196571   196389      -    no    1 0.33   5.3   12.8       1.3 ?   U2 spliceosomal RNA
U2                   RF00004   ENA|BK006936|BK006936.2 -          cm        1      193   681858   681747      -    no    1 0.33   0.1   98.7   1.2e-20 !   U2 spliceosomal RNA
5S_rRNA              RF00001   ENA|BK006937|BK006937.2 -          cm        1      119      761      644      -    no    1 0.41   0.3   14.1       2.4 ?   5S ribosomal RNA
U2                   RF00004   ENA|BK006937|BK006937.2 -          cm        1      193   229986   229885      -    no    1 0.32   0.1   11.1       4.7 ?   U2 spliceosomal RNA
#
# Program:         cmscan
# Version:         1.1.5 (Sep 2023)
# Pipeline mode:   SCAN
# Query file:      BK006935_6_7.fasta
# Target file:     IRES_5S_U2.cm
# Option settings: cmscan --tblout IRES_5S_U2_Yeast-cmscan.tbl IRES_5S_U2.cm BK006935_6_7.fasta 
# Current dir:     /analysis/BioPython/Testing
# Date:            Thu Sep 12 05:06:09 2024
# [ok]
```

## 1) `dnakit.alignment.biopython_p18_alignment_functions.read_alignment_records()`

- **作用：** 读取其他比对读取函数未覆盖的比对格式。
- **处理方法：** 对 `a2m`、`bigmaf`（bigMaf，二进制 bigBed 文件）、`bigpsl`（bigPsl）、`chain`（UCSC chain）、`emboss`（EMBOSS needle/water/matcher 输出）、`exonerate`（自动识别 text、cigar、vulgar 输出）、`mauve`（XMFA）、`msf`、`nexus`、`stockholm` 或 `tabular`（用 `-outfmt 7` 写出的 BLAST 或 FASTA 表格输出）调用 `Bio.Align.parse()`。块坐标为 0 起始、半开区间，每个序列一行。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.read_alignment_records(path, *, file_format='nexus', max_alignments=10000)`。
- **输入：** 本地文件路径、`file_format` 与 `max_alignments`（默认 10000）。
- **输出：** `alignment_count` 与 `alignments`；每个比对含 `sequence_ids`、`sequence_lengths`、`aligned_columns`、`score`、`coordinates`、文件保存序列时的带空位 `aligned_rows`、`annotations` 与 `column_annotations`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.biopython_p18_alignment_functions import (
    read_alignment_records,
    write_aligned_sequences,
)

path = Path(tempfile.mkdtemp()) / "demo.nex"
write_aligned_sequences({"s1": "ACGT-ACGTA", "s2": "AC-TTACGAA"}, path, file_format="nexus")
result = read_alignment_records(path, file_format="nexus")
alignment = result["alignments"][0]
print(result["alignment_count"], alignment["sequence_ids"], alignment["aligned_columns"])
print(alignment["coordinates"])
print(alignment["aligned_rows"])
```

- **示例输出：**

```text
1 ['s1', 's2'] 10
[[0, 2, 3, 4, 4, 9], [0, 2, 2, 3, 4, 9]]
{'s1': 'ACGT-ACGTA', 's2': 'AC-TTACGAA'}
```

## 2) `dnakit.alignment.biopython_p18_alignment_functions.write_aligned_sequences()`

- **作用：** 把一个多序列比对（序列名到带空位序列）写成新文件。
- **处理方法：** 由带空位的行构造 `Bio.Align.Alignment`（`Alignment.parse_printed_alignment()`），再用 `Bio.Align.write()` 写成 `a2m`、`clustal`、`fasta`、`maf`、`mauve`（XMFA）、`nexus`、`phylip` 或 `stockholm` 格式。`molecule_type`（`DNA` 或 `RNA`）供 NEXUS 写出器使用。PHYLIP 名称的填充与长度限制由 Biopython 写出器决定。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.write_aligned_sequences(aligned_sequences, output_path, *, file_format='clustal', molecule_type='DNA')`。
- **输入：** 名称到带空位序列的映射（长度相等、名称无空格）、新的输出路径、`file_format` 与 `molecule_type`。
- **输出：** `output_path`、`file_format`、`sequence_count`、`aligned_columns` 和 `alignment_count`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.biopython_p18_alignment_functions import write_aligned_sequences

out = Path(tempfile.mkdtemp()) / "demo.phy"
rows = {"s1": "ACGT-ACGTA", "s2": "AC-TTACGAA", "s3": "ACGTTACG-A"}
result = write_aligned_sequences(rows, out, file_format="phylip")
print(result["sequence_count"], result["aligned_columns"], result["alignment_count"])
print(out.read_text(encoding="utf-8"))
```

- **示例输出：**

```text
3 10 1
3 10
s1        ACGT-ACGTA
s2        AC-TTACGAA
s3        ACGTTACG-A
```

## 3) `dnakit.alignment.biopython_p18_alignment_functions.convert_alignment_file()`

- **作用：** 用 `Bio.Align` 把比对文件转换为 A2M、bigBed、bigMaf、bigPsl、chain、Clustal、Exonerate、MAF、Mauve、PHYLIP、PSL 或 Stockholm。
- **处理方法：** 用 `Bio.Align.parse()` 解析源文件（源格式：a2m、bigbed、bigmaf、bigpsl、chain、clustal、emboss、exonerate、fasta、maf、mauve、msf、nexus、phylip、psl、stockholm、tabular），再用 `Bio.Align.write()` 写出。Exonerate 目标由 `exonerate_variant` 选择 `vulgar` 或 `cigar`。会打印序列的文本目标（a2m、clustal、maf、mauve、phylip、stockholm）要求源文件保存序列；坐标型目标（bigBed、bigMaf、bigPsl、chain、Exonerate、PSL）要求 Biopython 能以该格式写出的成对比对源。Biopython 写出器从其他格式的比对中读取的属性（目标大小、注释、列状态、Mauve 标识符）由函数补全，bigBed 系列还会按目标名和起点排序。Biopython 拒绝某种组合时，函数带着其消息报错并删除不完整文件。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.convert_alignment_file(source, output_path, *, source_format='psl', target_format='chain', exonerate_variant='vulgar', max_alignments=10000)`。
- **输入：** 源路径、新的输出路径、`source_format`、`target_format`、`exonerate_variant` 与 `max_alignments`。
- **输出：** `alignment_count`、`bytes_written`、`source_format` 与 `target_format`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dnakit.alignment.biopython_p18_alignment_functions import convert_alignment_file

aligner = Align.PairwiseAligner(
    mode="global", match_score=2, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1
)
target = SeqRecord(Seq("ACGTACGTTGCA"), id="chr1")
query = SeqRecord(Seq("ACGTCGTTTGCA"), id="read1")
work = Path(tempfile.mkdtemp())
Align.write(aligner.align(target, query)[0], work / "demo.psl", "psl")
result = convert_alignment_file(work / "demo.psl", work / "demo.chain", source_format="psl", target_format="chain")
print(result["alignment_count"], result["target_format"])
print((work / "demo.chain").read_text(encoding="utf-8"))
```

- **示例输出：**

```text
1 chain
chain 0 chr1 12 + 0 12 read1 12 + 0 12
4	1	0
2	0	1
5
```

## 4) `dnakit.alignment.biopython_p18_alignment_functions.search_maf_index()`

- **作用：** 找出 MAF 文件中与参考序列区间重叠的比对块。
- **处理方法：** 在 `target_seqname`（例如 `mm9.chr10`）上构建 `Bio.AlignIO.MafIO.MafIndex` 并调用 `search()`。索引保存在 SQLite 文件 `index_path` 中，已存在时直接复用；不给 `index_path` 则使用临时索引并在结束时删除。`starts` 为 0 起始，`ends` 为半开区间，遵循 UCSC 约定。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.search_maf_index(maf_path, target_seqname, starts, ends, *, index_path=None, max_blocks=1000)`。
- **输入：** MAF 路径、`target_seqname`、`starts`、`ends`，以及可选的 `index_path` 与 `max_blocks`。
- **输出：** `indexed_blocks`、`block_count` 与 `blocks`；每个块列出 `rows`（`id`、带空位的 `sequence`、MAF 的 `annotations`，如 start、size、strand、srcSize）。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import search_maf_index

result = search_maf_index("demo.maf", "hg.chr1", [102], [108])
print(result["indexed_blocks"], result["block_count"])
for row in result["blocks"][0]["rows"]:
    print(row["id"], row["sequence"], row["annotations"]["start"])
```

- **示例输出：**

```text
2 1
hg.chr1 ACGTACGTAC 100
mm.chr2 ACGTA-GTAC 200
```

## 5) `dnakit.alignment.biopython_p18_alignment_functions.splice_maf_index()`

- **作用：** 提取若干参考序列片段（例如外显子）的拼接比对。
- **处理方法：** 在与 `search_maf_index()` 相同的索引上调用 `MafIndex.get_spliced()`：把 `target_seqname` 上 `starts`/`ends`（0 起始、半开）指定的片段从 MAF 比对中切出并拼接，得到每个物种一行的多序列比对。`strand` 为 1 或 -1；-1 时对结果取反向互补。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.splice_maf_index(maf_path, target_seqname, starts, ends, *, strand=1, index_path=None)`。
- **输入：** MAF 路径、`target_seqname`、`starts`、`ends`、`strand`，以及可选的 `index_path`。
- **输出：** `alignment_length`、`segment_count`、`strand` 与 `rows`（序列名到拼接后带空位序列）。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import splice_maf_index

forward = splice_maf_index("demo.maf", "hg.chr1", [100, 120], [104, 124])
print(forward["alignment_length"], forward["rows"])
reverse = splice_maf_index("demo.maf", "hg.chr1", [100, 120], [104, 124], strand=-1)
print(reverse["rows"])
```

- **示例输出：**

```text
8 {'hg.chr1': 'ACGTGGTT', 'mm.chr2': 'ACGTGGTT'}
{'hg.chr1': 'AACCACGT', 'mm.chr2': 'AACCACGT'}
```

## 6) `dnakit.alignment.biopython_p18_alignment_functions.read_search_results()`

- **作用：** 用 `Bio.SearchIO` 读取检索程序报告（HMMER、Infernal 及其他格式）。
- **处理方法：** 对 `hmmer2-text`、`hmmer3-text`、`hmmer3-tab`（`--tblout`）、`hmmscan3-domtab`、`hmmsearch3-domtab`、`phmmer3-domtab`（`--domtblout`）、`infernal-tab`、`infernal-text`、`blast-tab`、`blast-xml`、`blat-psl`、`exonerate-text`、`exonerate-cigar`、`exonerate-vulgar` 或 `fasta-m10` 调用 `Bio.SearchIO.parse()`。`comments=True` 读取用 `-outfmt 7` 写出的 BLAST 表格文件。HSP 范围为 0 起始、半开区间，与 Bio.SearchIO 的存储方式相同。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.read_search_results(path, *, file_format='hmmer3-text', comments=False, max_queries=10000)`。
- **输入：** 本地报告路径、`file_format`、`comments` 与 `max_queries`。
- **输出：** `query_count`、`hit_count` 与 `queries`；每个查询含程序、目标、版本、`hits`（id、描述，以及格式提供时的 E 值与比特分数）和 `hsps`（范围、分数、链方向及该格式特有的属性）。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import read_search_results

result = read_search_results("hits.domtab", file_format="hmmscan3-domtab")
print(result["query_count"], result["hit_count"])
for query in result["queries"]:
    for hit in query["hits"]:
        first = hit["hsps"][0]
        print(hit["id"], hit["hsp_count"], first["evalue"], first["query_range"])
```

- **示例输出：**

```text
2 4
Globin 1 1.6e-22 [6, 112]
Xpo1 2 1.2e-33 [109, 271]
IBN_N 2 0.036 [35, 98]
Rac1 2 1.3 [88, 169]
```

## 7) `dnakit.alignment.biopython_p18_alignment_functions.index_search_results()`

- **作用：** 为检索输出文件建立按查询 id 随机访问的索引，并返回选定的查询。
- **处理方法：** 不给 `index_path` 时建立内存偏移索引（`Bio.SearchIO.index()`）；给出时索引保存在 SQLite 文件中（`Bio.SearchIO.index_db()`），不存在则创建，存在则复用。通过索引读取 `query_ids` 指定的查询。Biopython 1.88 无法为 HMMER 表格格式（`hmmer3-tab` 与三种 `domtab` 格式）建立索引，这些格式会被明确拒绝。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.index_search_results(path, *, file_format='hmmer3-text', index_path=None, query_ids=None, comments=False)`。
- **输入：** 报告路径、`file_format`，以及可选的 `index_path`、`query_ids` 与 `comments`。
- **输出：** `query_count`、排序后的 `query_ids`、`index_path` 以及完整的所请求 `queries`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.alignment.biopython_p18_alignment_functions import index_search_results

database = Path(tempfile.mkdtemp()) / "hits.sqlite"
result = index_search_results(
    "cm.tbl", file_format="infernal-tab", index_path=database,
    query_ids=["ENA|BK006936|BK006936.2"],
)
print(result["query_count"], result["query_ids"])
print(result["queries"][0]["id"], result["queries"][0]["hit_count"], database.exists())
```

- **示例输出：**

```text
3 ['ENA|BK006935|BK006935.2', 'ENA|BK006936|BK006936.2', 'ENA|BK006937|BK006937.2']
ENA|BK006936|BK006936.2 1 True
```

## 8) `dnakit.alignment.biopython_p18_alignment_functions.filter_search_results()`

- **作用：** 过滤检索报告中的命中和 HSP。
- **处理方法：** 先用 `QueryResult.hsp_filter()` 保留 E 值不超过 `max_evalue`、比特分数不低于 `min_bitscore` 且比对跨度不小于 `min_alignment_span` 的 HSP；失去全部 HSP 的命中随之消失。再用 `QueryResult.hit_filter()` 保留 HSP 数不少于 `min_hsps_per_hit` 的命中，给出 `hit_ids` 时还要求 id 在该列表中。未设置的条件被忽略，至少要给一个。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.filter_search_results(path, *, file_format='hmmer3-text', max_evalue=None, min_bitscore=None, min_alignment_span=None, min_hsps_per_hit=None, hit_ids=None, comments=False)`。
- **输入：** 报告路径、`file_format` 与上述条件。
- **输出：** `hits_before`、`hits_after` 与剩余的 `queries`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import filter_search_results

strict = filter_search_results("hits.domtab", file_format="hmmscan3-domtab", max_evalue=1e-5)
print(strict["hits_before"], strict["hits_after"])
print([hit["id"] for query in strict["queries"] for hit in query["hits"]])
multi = filter_search_results("hits.domtab", file_format="hmmscan3-domtab", min_hsps_per_hit=2)
print(multi["hits_after"])
```

- **示例输出：**

```text
4 2
['Globin', 'Xpo1']
3
```

## 9) `dnakit.alignment.biopython_p18_alignment_functions.sort_search_results()`

- **作用：** 对检索报告中每个查询的命中以及每个命中的 HSP 排序。
- **处理方法：** 用 `QueryResult.sort()` 与 `Hit.sort()` 按命名键排序：`evalue`（升序，无 E 值的命中排最后）、`bitscore`、`id`、`length`（格式保存时的命中序列长度）或 `hsp_count`。`reverse` 反转顺序。`sort_hsps` 为真且键适用于 HSP（evalue、bitscore）时，HSP 按同一键排序。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.sort_search_results(path, *, file_format='hmmer3-text', sort_by='evalue', reverse=False, sort_hsps=True, comments=False)`。
- **输入：** 报告路径、`file_format`、`sort_by`、`reverse` 与 `sort_hsps`。
- **输出：** `sort_by`、`reverse` 与排序后的 `queries`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import sort_search_results

result = sort_search_results("hits.domtab", file_format="hmmscan3-domtab", sort_by="bitscore", reverse=True)
query = result["queries"][1]
print([(hit["id"], hit["bitscore"]) for hit in query["hits"]])
print([hsp["bitscore"] for hsp in query["hits"][0]["hsps"]])
```

- **示例输出：**

```text
[('Xpo1', 116.6), ('IBN_N', 16.9), ('Rac1', 11.7)]
[116.1, -1.8]
```

## 10) `dnakit.alignment.biopython_p18_alignment_functions.fogsaa_global_alignment()`

- **作用：** 用 FOGSAA 算法对两条序列做全局比对。
- **处理方法：** 使用 `Bio.Align.PairwiseAligner(mode='fogsaa')`，即针对最优全局比对的分支定界搜索。分数默认取 Biopython 的值（匹配 1、错配 0、无空位罚分）；要得到有意义的 DNA 比对，请传入明确的分数。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.fogsaa_global_alignment(sequence_a, sequence_b, *, match_score=1.0, mismatch_score=0.0, open_gap_score=0.0, extend_gap_score=0.0, max_alignments=1)`。
- **输入：** `sequence_a`、`sequence_b`、`match_score`、`mismatch_score`、`open_gap_score`、`extend_gap_score` 与 `max_alignments`。
- **输出：** `score`、`alignment_count`、`parameters` 与 `alignments`（名为 `sequence_a`、`sequence_b` 的带空位行以及块坐标）。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import fogsaa_global_alignment

result = fogsaa_global_alignment(
    "ACGTACGTTGCA", "ACGTCGTTTGCA",
    match_score=2, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1,
)
print(result["score"], result["alignment_count"])
rows = result["alignments"][0]["aligned_rows"]
print(rows["sequence_a"])
print(rows["sequence_b"])
```

- **示例输出：**

```text
18.0 1
ACGTACGTT-GCA
ACGT-CGTTTGCA
```

## 11) `dnakit.alignment.biopython_p18_alignment_functions.align_nucleotide_to_protein()`

- **作用：** 用密码子感知比对器把编码核苷酸序列与其蛋白质比对。
- **处理方法：** `Bio.Align.CodonAligner` 把每个氨基酸与其密码子配对，并容忍核苷酸序列中的移码和内含子；`anchor_length` 是在空位之后寻找下一个密码子所用的种子长度。接受 DNA 或 RNA（RNA 按 DNA 处理）。`codon_table` 是 NCBI 遗传密码表编号。分数为比对上的密码子数。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.align_nucleotide_to_protein(protein, nucleotides, *, codon_table=1, anchor_length=10, max_alignments=10)`。
- **输入：** `protein`、`nucleotides`、`codon_table`、`anchor_length` 与 `max_alignments`。
- **输出：** `score`、`alignment_count` 与带块坐标的 `alignments`（蛋白以残基计，核苷酸以碱基计）。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import align_nucleotide_to_protein

plain = align_nucleotide_to_protein("MAIVMG", "ATGGCCATTGTAATGGGC")
print(plain["score"], plain["alignments"][0]["coordinates"])
shifted = align_nucleotide_to_protein(
    "MLGSNESRVCPKLSPS", "ATGCTGGGCTCGAACGAGTCCGTGTATGCCCTAAGCTGAGCCCGTCG"
)
print(shifted["score"], shifted["alignments"][0]["coordinates"])
```

- **示例输出：**

```text
6.0 [[0, 6], [0, 18]]
13.0 [[0, 7, 7, 16], [0, 21, 20, 47]]
```

## 12) `dnakit.alignment.biopython_p18_alignment_functions.sort_alignment_rows()`

- **作用：** 对比对的各行排序。
- **处理方法：** 构建 `Alignment` 并用命名键调用 `Alignment.sort()`：`id`（Biopython 默认）、`sequence`（去空位后的字母）、`gc_fraction`（去空位序列的 GC 比例）、`ungapped_length` 或 `gap_count`。并列时保持输入顺序。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.sort_alignment_rows(aligned_sequences, *, sort_by='id', reverse=False)`。
- **输入：** 名称到带空位序列的映射、`sort_by` 与 `reverse`。
- **输出：** `order`（新顺序下的名称）以及按该顺序排列的 `aligned_sequences`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import sort_alignment_rows

rows = {"r1": "ACGT-ACGTA", "r2": "AC-TTACGAA", "r3": "ACGTTACG-A"}
print(sort_alignment_rows(rows, sort_by="id", reverse=True)["order"])
print(sort_alignment_rows(rows, sort_by="gc_fraction", reverse=True)["order"])
print(sort_alignment_rows(rows, sort_by="gap_count", reverse=True)["order"])
```

- **示例输出：**

```text
['r3', 'r2', 'r1']
['r1', 'r3', 'r2']
['r1', 'r2', 'r3']
```

## 13) `dnakit.alignment.biopython_p18_alignment_functions.alignment_substitution_counts()`

- **作用：** 统计比对各行之间的字母替换数。
- **处理方法：** 读取 `Alignment.substitutions`：(行字母, 列字母) 项是第一个字母与第二个字母对齐的列数（多行比对时按 Biopython 的方式对称化）。空位列不计入。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.alignment_substitution_counts(aligned_sequences)`。
- **输入：** 名称到带空位序列的映射（至少两行）。
- **输出：** `alphabet` 与方阵 `matrix`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_substitution_counts

result = alignment_substitution_counts({"a": "ACGTTACGA", "b": "ACGATACGA", "c": "ACCTTACGA"})
print(result["alphabet"])
for letter, row in zip(result["alphabet"], result["matrix"]):
    print(letter, row)
```

- **示例输出：**

```text
ACGT
A [9.0, 0.0, 0.0, 1.0]
C [0.0, 6.0, 0.0, 0.0]
G [0.0, 2.0, 4.0, 0.0]
T [1.0, 0.0, 0.0, 4.0]
```

## 14) `dnakit.alignment.biopython_p18_alignment_functions.alignment_column_frequencies()`

- **作用：** 统计比对每一列中各字母的个数。
- **处理方法：** 读取 `Alignment.frequencies`：对每个字母（空位用 `-`），给出每一列中带有该字母的行数。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.alignment_column_frequencies(aligned_sequences)`。
- **输入：** 名称到带空位序列的映射（至少两行）。
- **输出：** `column_count`、`row_count` 与 `frequencies`（字母到逐列计数列表）。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_column_frequencies

rows = {"r1": "ACGT-ACGTA", "r2": "AC-TTACGAA", "r3": "ACGTTACG-A"}
result = alignment_column_frequencies(rows)
print(result["row_count"], result["column_count"])
for letter, counts in result["frequencies"].items():
    print(letter, counts)
```

- **示例输出：**

```text
3 10
- [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
A [3.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 1.0, 3.0]
C [0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0]
G [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0]
T [0.0, 0.0, 0.0, 3.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0]
```

## 15) `dnakit.alignment.biopython_p18_alignment_functions.alignment_identity_counts()`

- **作用：** 统计比对的一致位点、错配和空位数。
- **处理方法：** 调用 `Alignment.counts()`。返回对齐列数、一致位点数、错配数、空位数，以及其按左、内部、右空位和按插入、缺失的细分。`wildcard`（单个字符，如 `N`）使该字母不参与计数；`scoring_matrix`（内置替换矩阵名，如 `NUC.4.4`）会增加 `substitution_score` 与 `positives` 数。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.alignment_identity_counts(aligned_sequences, *, wildcard=None, scoring_matrix=None)`。
- **输入：** 名称到带空位序列的映射，以及 `wildcard` 或 `scoring_matrix`。
- **输出：** `row_count`、`column_count` 与 `counts`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_identity_counts

pair = {"a": "ACGT-ACGTA", "b": "AC-TTACGAA"}
counts = alignment_identity_counts(pair)["counts"]
print(counts["aligned"], counts["identities"], counts["mismatches"], counts["gaps"], counts["internal_gaps"])
scored = alignment_identity_counts(pair, scoring_matrix="NUC.4.4")["counts"]
print(scored["positives"], scored["substitution_score"])
```

- **示例输出：**

```text
8 7 1 2 2
7 31.0
```

## 16) `dnakit.alignment.biopython_p18_alignment_functions.alignment_coordinate_maps()`

- **作用：** 在比对列与序列位置之间建立映射。
- **处理方法：** 覆盖 `Alignment.aligned`、`Alignment.indices` 与 `Alignment.inverse_indices`。`indices` 给出每一行每一列在去空位序列中的 0 起始位置，空位处为 -1；`inverse_indices` 给出每个去空位位置所在的列；`aligned_blocks`（仅两行时）按序列列出对齐块区间，多于两行时为 null。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.alignment_coordinate_maps(aligned_sequences)`。
- **输入：** 名称到带空位序列的映射（至少两行）。
- **输出：** `indices`、`inverse_indices`、`aligned_blocks` 与 `column_count`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import alignment_coordinate_maps

result = alignment_coordinate_maps({"a": "AC-GT", "b": "ACAG-"})
print(result["indices"])
print(result["inverse_indices"])
print(result["aligned_blocks"])
```

- **示例输出：**

```text
{'a': [0, 1, -1, 2, 3], 'b': [0, 1, 2, 3, -1]}
{'a': [0, 1, 3, 4], 'b': [0, 1, 2, 3]}
[[[0, 2], [2, 3]], [[0, 2], [3, 4]]]
```

## 17) `dnakit.alignment.biopython_p18_alignment_functions.parse_printed_alignment()`

- **作用：** 由打印出的比对行推断无空位序列和块坐标。
- **处理方法：** 调用 `Alignment.parse_printed_alignment()`：每个序列一行带空位文本，所有行列数相同，空位用 `-`。返回去空位序列和坐标矩阵（行是序列、列是块边界，0 起始），二者共同定义一个 `Alignment`。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.parse_printed_alignment(aligned_lines)`。
- **输入：** 带空位行的列表。
- **输出：** `sequences` 与 `coordinates`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import parse_printed_alignment

result = parse_printed_alignment(["TAGGCATACGTG", "AACG--TACGT-", "-ACGCATACTTG"])
print(result["sequences"])
for row in result["coordinates"]:
    print(row)
```

- **示例输出：**

```text
['TAGGCATACGTG', 'AACGTACGT', 'ACGCATACTTG']
[0, 1, 4, 6, 11, 12]
[0, 1, 4, 4, 9, 9]
[0, 0, 3, 5, 10, 11]
```

## 18) `dnakit.alignment.biopython_p18_alignment_functions.merge_alignments_by_reference()`

- **作用：** 把共享同一参考序列的成对比对合并成多序列比对。
- **处理方法：** 调用 `Alignment.from_alignments_with_same_reference()`。每个查询对应 `{"target": 带空位的参考行, "query": 带空位的查询行}`；去掉空位后所有参考序列必须相同。结果中参考序列排在最前，查询插入字母处会增加列。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.merge_alignments_by_reference(reference_alignments)`。
- **输入：** 查询名到 `target`/`query` 带空位行对的映射。
- **输出：** `row_count`、`column_count` 与 `aligned_sequences`（参考序列名为 `reference`）。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import merge_alignments_by_reference

result = merge_alignments_by_reference(
    {
        "q1": {"target": "ACGT", "query": "AC-T"},
        "q2": {"target": "ACGT-", "query": "A-GTT"},
    }
)
print(result["row_count"], result["column_count"])
for name, row in result["aligned_sequences"].items():
    print(name, row)
```

- **示例输出：**

```text
3 5
reference ACGT-
q1 AC-T-
q2 A-GTT
```

## 19) `dnakit.alignment.biopython_p18_alignment_functions.map_alignment_through()`

- **作用：** 让一个比对通过其查询序列的一个或多个比对进行映射（坐标转换）。
- **处理方法：** `alignment` 为 `{"target": 带空位行, "query": 带空位行}`，例如染色体对转录本。`through` 列出同样形式的比对，其 target 必须是第一个比对的 query（去空位后），例如比对到转录本上的 reads。单个比对用 `Alignment.map()`，多个比对用 `Alignment.mapall()`。结果是映射后的序列对第一个 target 的比对。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.map_alignment_through(alignment, through)`。
- **输入：** 第一个比对以及要经过的比对列表。
- **输出：** `method`（`map` 或 `mapall`）、`coordinates` 与带空位的 `aligned_rows`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import map_alignment_through

chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
transcript_alignment = {"target": chromosome, "query": "--------CCCCCCC-----------GGGGGG--------"}
read_alignment = {"target": "CCCCCCCGGGGGG", "query": "---CCCCGGGG--"}
result = map_alignment_through(transcript_alignment, [read_alignment])
print(result["method"], result["coordinates"])
print(result["aligned_rows"]["target"])
print(result["aligned_rows"]["query"])
```

- **示例输出：**

```text
map [[11, 15, 26, 30], [0, 4, 4, 8]]
CCCCAAAAAAAAAAAGGGG
CCCC-----------GGGG
```

## 20) `dnakit.alignment.biopython_p18_alignment_functions.find_alignment_columns()`

- **作用：** 列出 DNA 比对中的恒定列或全空位列。
- **处理方法：** `kind="constant"` 调用 `Nexus.constant()`：所有（未删除）序列一致的列，`?` 视为相容，简并码视为碱基集合。`kind="gap_only"` 调用 `Nexus.gaponly()`：只含空位的列，`include_missing` 为真时 `?` 也算空位。`delete_taxa` 排除序列，`exclude_sites` 在恒定列搜索中排除列。列号从 0 开始。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.find_alignment_columns(aligned_sequences, *, kind='constant', delete_taxa=None, exclude_sites=None, include_missing=False)`。
- **输入：** 名称到带空位 DNA 的映射、`kind`、`delete_taxa`、`exclude_sites` 与 `include_missing`。
- **输出：** `kind`、`columns`、`found` 与 `column_count`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import find_alignment_columns

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
print(find_alignment_columns(dna)["columns"])
print(find_alignment_columns(dna, kind="gap_only")["columns"])
print(find_alignment_columns(dna, delete_taxa=["t3"])["columns"])
```

- **示例输出：**

```text
[0, 1, 2, 4, 5, 7]
[4]
[0, 1, 2, 3, 4, 5, 7]
```

## 21) `dnakit.alignment.biopython_p18_alignment_functions.site_character_states()`

- **作用：** 汇总比对各列的字符状态。
- **处理方法：** 对 `sites` 中的每一列（为 null 时取全部列，0 起始）调用 `Nexus.cstatus()`：返回未删除序列中出现的碱基的排序列表，简并码会展开。`narrow=True` 遵循 PAUP（`A C ?` 得到 `A C`；全是 `?` 则得到 `?`）；`narrow=False` 把 `?` 和空位展开为所有状态。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.site_character_states(aligned_sequences, *, sites=None, delete_taxa=None, narrow=True)`。
- **输入：** 名称到带空位 DNA 的映射、`sites`、`delete_taxa` 与 `narrow`。
- **输出：** `sites`（`{site, states}` 列表）与 `column_count`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import site_character_states

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
for item in site_character_states(dna, sites=[3, 4, 6])["sites"]:
    print(item["site"], item["states"])
print(site_character_states(dna, sites=[4], narrow=False)["sites"][0]["states"])
```

- **示例输出：**

```text
3 ['A', 'T']
4 ['?']
6 ['C', 'T']
['-', 'A', 'C', 'G', 'T']
```

## 22) `dnakit.alignment.biopython_p18_alignment_functions.subset_alignment()`

- **作用：** 从 DNA 比对中删除序列并排除列。
- **处理方法：** 调用 `Nexus.crop_matrix()`。`delete_taxa` 指定要删除的序列，`exclude_sites` 指定要移除的 0 起始列。两者至少给一个，且至少保留一条序列。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.subset_alignment(aligned_sequences, *, delete_taxa=None, exclude_sites=None)`。
- **输入：** 名称到带空位 DNA 的映射、`delete_taxa` 与 `exclude_sites`。
- **输出：** `row_count`、`column_count` 与 `aligned_sequences`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import subset_alignment

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
result = subset_alignment(dna, delete_taxa=["t2"], exclude_sites=[0, 1])
print(result["row_count"], result["column_count"])
print(result["aligned_sequences"])
```

- **示例输出：**

```text
2 6
{'t1': 'GT-ACG', 't3': 'GA-ACG'}
```

## 23) `dnakit.alignment.biopython_p18_alignment_functions.insert_alignment_gaps()`

- **作用：** 向 DNA 比对中插入空位列。
- **处理方法：** 调用 `Nexus.insert_gap()`：在 0 起始的列 `position` 之前，给每条序列加入 `count` 列 `-`。插入点位于一段空位内部时，新空位放在该段右端，`left_greedy=True` 时放在左端。
- **API：** `dnakit.alignment.biopython_p18_alignment_functions.insert_alignment_gaps(aligned_sequences, position, *, count=1, left_greedy=False)`。
- **输入：** 名称到带空位 DNA 的映射、`position`、`count` 与 `left_greedy`。
- **输出：** `column_count` 与 `aligned_sequences`。
- **示例代码：**

```python
from dnakit.alignment.biopython_p18_alignment_functions import insert_alignment_gaps

dna = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}
result = insert_alignment_gaps(dna, 2, count=2)
print(result["column_count"])
for name, row in result["aligned_sequences"].items():
    print(name, row)
```

- **示例输出：**

```text
10
t1 AC--GT-ACG
t2 AC--GT-ATG
t3 AC--GA-ACG
```
