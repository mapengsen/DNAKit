# G18o 功能：Biopython、scikit-bio、cutadapt 与 umi_tools 的评审包装

本页记录评审组 G18o 新增的 26 个包装函数，覆盖 `sequence` 与 `sequencing_qc` 两个区域。这些函数都从对应的区域包导出，例如 `from dnakit.sequence import write_motifs`。

- 离线运行：所有函数只读取调用参数或本地文件（SAM、BAM、FASTA、QUAL、BLAST 表格与 Newick 字符串），不访问网络。
- 写文件的函数只创建新文件；输出路径已存在时抛出 `ConfigurationError`，因此示例反复运行也不会冲突。
- 示例所需的输入以 `示例文件` 标签给出，示例代码在运行目录中读取它们；需要中间文件的示例使用临时目录。
- pyslow5 1.5.0 遇到损坏的 SLOW5/BLOW5 文件时会使解释器崩溃。本组没有为 pyslow5 新增函数；现有 `read_slow5_signal` 等后端在打开文件前检查 SLOW5/BLOW5 魔数字节，调用时应保留这一检查。
- parasail 的 `ssw_init` 在 score_size 为 16 时会段错误。本组没有包装它；将来包装时必须先对 score_size 做保护。

## 台账状态（88 行 other 与 cap）

| 包 | 原有 other/cap | 包装 | 其他状态 |
| --- | ---: | ---: | --- |
| Bio | 39 | 13 行（14 个函数，含 `write_gaf_file`） | family-variant 11、duplicate-of 1、not-dna 6、internal 6、deprecated 2 |
| skbio | 10 | 5 | internal 4、not-dna 1 |
| parasail | 13 | 0 | internal 7、not-dna 6（CPU 探测与计时器） |
| mappy | 1 | 0 | internal 1（`mappy.verbose`） |
| pyhmmer | 1 | 0 | internal 1（`hmmpress`） |
| cutadapt | 11 | 1 | cli 4、internal 5、duplicate-of 1 |
| umi_tools | 11 | 6 | family-variant 1、cli 1、internal 3 |
| pyfastx | 1 | 0 | internal 1（`version`） |
| pyslow5 | 1 | 0 | already-wrapped 1（现有 `read_slow5_signal` 等） |

未包装的主要原因：`not-dna` 为 Entrez、MEDLINE、SCAD 与 SFF 元数据读写，以及 parasail 的 CPU 探测；`internal` 为内部辅助、全局设置函数（会改变进程状态）与报告格式化；`cli` 为命令行报告与 RSEM 脚本的步骤；`deprecated` 为已退役的 UniGene。

## 1) `write_pairwise_alignment_file`

- **作用：** 用 Biopython 对齐两条 DNA 序列，并把最佳比对写入一个新文件。
- **处理方法：** 打分固定为匹配 2、错配 -1、gap open -2、gap extend -1，默认全局比对；`mode="local"` 为局部比对。格式可为 fasta、clustal、maf、phylip、psl、sam 或 bed。输出文件必须不存在。
- **API：** `write_pairwise_alignment_file(query[必须], target[必须], output_path[必须], file_format[可选], mode[可选], query_id[可选], target_id[可选], verbose[可选])`
- **输入：** 两条 DNA 字符串（允许 IUPAC 符号）、新文件路径、格式名、比对模式与序列名。
- **输出：** 字典，含 `score`、`aligned_query`、`aligned_target`、`alignment_count`、`output_path` 与 `provenance`。

**示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import write_pairwise_alignment_file

out = Path(tempfile.mkdtemp()) / "pair.fa"
result = write_pairwise_alignment_file("ACGTACGT", "ACGAACGT", out, verbose=False)
print(result["score"], result["aligned_query"], result["aligned_target"])
print(out.read_text().strip())
```

**示例输出：**

```text
13.0 ACGTACGT ACGAACGT
>query
ACGTACGT
>target
ACGAACGT
```

## 2) `convert_search_results`

- **作用：** 把本地 BLAST、BLAT 或 HMMER 搜索结果转换为另一种格式。
- **处理方法：** 调用 Bio.SearchIO.convert。输入格式与输出格式必须在白名单中；上游不能从某种输入写出某种输出时（例如 blast-tab 无法写出 blast-xml），返回 `ConfigurationError`。
- **API：** `convert_search_results(source[必须], source_format[必须], output_path[必须], output_format[必须], verbose[可选])`
- **输入：** 搜索结果文件（例如 BLAST 表格，无表头）、输入格式名与新输出文件路径。
- **输出：** 字典，含 `query_count`、`hit_count`、`hsp_count`、`fragment_count`、`output_path` 与 `provenance`。

**示例文件：** `hits.tsv`

```text
q1	t1	100.00	8	0	0	1	8	1	8	1e-05	16.0
q2	t2	87.50	8	1	0	1	8	1	8	2e-03	12.5
```

**示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import convert_search_results

out = Path(tempfile.mkdtemp()) / "converted.tsv"
result = convert_search_results("hits.tsv", "blast-tab", out, "blast-tab", verbose=False)
print(result["query_count"], result["hit_count"], result["hsp_count"])
print(out.read_text().strip())
```

**示例输出：**

```text
2 2 2
q1	t1	100.00	8	0	0	1	8	1	8	1e-05	16.0
q2	t2	87.50	8	1	0	1	8	1	8	0.002	12.5
```

## 3) `codon_back_table`

- **作用：** 返回遗传密码表中每个氨基酸对应的一个密码子，终止子记为 `*`。
- **处理方法：** 调用 Bio.Data.CodonTable.make_back_table。该表只保留每个氨基酸排序后的最后一个密码子，因此结果是确定的，但不代表全部同义密码子。
- **API：** `codon_back_table(genetic_code[可选], stop_codon[可选], verbose[可选])`
- **输入：** NCBI 遗传密码表编号（1 到 33 中的无歧义 DNA 表）与一个属于该表的终止密码子。
- **输出：** 字典，含 `table_name`、`codon_count`、`back_table`（氨基酸到密码子的映射）与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import codon_back_table

result = codon_back_table(1, verbose=False)
back = result["back_table"]
print(result["table_name"], result["codon_count"])
print(back["M"], back["L"], back["*"])
```

**示例输出：**

```text
Standard 61
ATG TTG TAA
```

## 4) `adam_consensus_tree`

- **作用：** 计算多棵 Newick 树的 Adam 一致树。
- **处理方法：** 每棵树先解析为 Bio.Phylo 树，再调用 Bio.Phylo.Consensus.adam_consensus。输出为上游结果，一致树分支长度由上游设为 0；拓扑冲突时以上游算法为准。
- **API：** `adam_consensus_tree(newick_trees[必须], verbose[可选])`
- **输入：** 一个由 Newick 字符串组成的序列，1 到 1000 棵树，且树的叶子名称一致。
- **输出：** 字典，含 `tree_count`、`newick` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import adam_consensus_tree

newick = "((A:1,B:1):1,(C:1,D:1):1);"
result = adam_consensus_tree([newick, newick], verbose=False)
print(result["tree_count"], result["newick"])
```

**示例输出：**

```text
2 ((C:1,D:1):0,(A:1,B:1):0):0;
```

## 5) `bootstrap_branch_support`

- **作用：** 在目标树上标注每个内部分支被多少棵重复树支持。
- **处理方法：** 调用 Bio.Phylo.Consensus.get_support。支持度为包含该分支的重复树所占百分比，写在 Newick 的内部节点上。例如 3 棵重复树中有 2 棵含某分支，则为 66.67。
- **API：** `bootstrap_branch_support(target_newick[必须], replicate_newicks[必须], verbose[可选])`
- **输入：** 目标树的 Newick 字符串与 1 到 1000 棵重复树，叶子名称须与目标树相同。
- **输出：** 字典，含 `replicate_count`、`scale`（固定为 percent）、`newick` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import bootstrap_branch_support

target = "((A:0.1,B:0.2):0.3,(C:0.1,D:0.2):0.3);"
replicates = [target, "((A:0.1,C:0.2):0.3,(B:0.1,D:0.2):0.3);", target]
result = bootstrap_branch_support(target, replicates, verbose=False)
print(result["replicate_count"], result["scale"])
print(result["newick"])
```

**示例输出：**

```text
3 percent
((A:0.1,B:0.2)66.67:0.3,(C:0.1,D:0.2)66.67:0.3)100.00:0;
```

## 6) `tree_to_graph_edges`

- **作用：** 把 Newick 树的节点与父到子的边以普通列表返回，无需 networkx 或 igraph。
- **处理方法：** 按前序遍历给节点编号（根为 0），边为每个父节点到其子节点。与 `Bio.Phylo.to_networkx` 和 `Bio.Phylo.to_igraph` 得到的边集合相同。
- **API：** `tree_to_graph_edges(newick[必须], verbose[可选])`
- **输入：** 一个 Newick 字符串。
- **输出：** 字典，含 `node_count`、`edge_count`、`tip_count`、`nodes`（编号、名称、分支长度、是否叶子）、`edges`（源、目标、分支长度）与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import tree_to_graph_edges

result = tree_to_graph_edges("((A:0.1,B:0.2):0.3,C:0.5);", verbose=False)
print(result["node_count"], result["edge_count"], result["tip_count"])
for edge in result["edges"]:
    print(edge["source"], edge["target"], edge["branch_length"])
```

**示例输出：**

```text
5 4 3
0 1 0.3
0 4 0.5
1 2 0.1
1 3 0.2
```

## 7) `complement_rna_sequence`

- **作用：** 返回 RNA 序列的互补链，保持 U 字母。
- **处理方法：** 调用 Bio.Seq.complement_rna（A 与 U、C 与 G 互补）。允许 IUPAC 符号；输入中出现 T 会被拒绝。
- **API：** `complement_rna_sequence(sequence[必须], verbose[可选])`
- **输入：** 一个 RNA 字符串（ACGU 与 IUPAC 符号）。
- **输出：** 字典，含 `input`、`sequence`、`length` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import complement_rna_sequence

print(complement_rna_sequence("ACGU", verbose=False)["sequence"])
```

**示例输出：**

```text
UGCA
```

## 8) `reverse_complement_rna_sequence`

- **作用：** 返回 RNA 序列的反向互补链，保持 U 字母。
- **处理方法：** 调用 Bio.Seq.reverse_complement_rna：先反向，再互补。
- **API：** `reverse_complement_rna_sequence(sequence[必须], verbose[可选])`
- **输入：** 一个 RNA 字符串（ACGU 与 IUPAC 符号）。
- **输出：** 字典，含 `input`、`sequence`、`length` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import reverse_complement_rna_sequence

print(reverse_complement_rna_sequence("AAGC", verbose=False)["sequence"])
```

**示例输出：**

```text
GCUU
```

## 9) `read_fasta_qual_pairs`

- **作用：** 读取成对的 FASTA 与 QUAL 文件，得到带 PHRED 质量值的序列。
- **处理方法：** 调用 Bio.SeqIO.QualityIO.PairedFastaQualIterator，按位置配对两个文件。质量值为整数。记录数超过 `max_records` 时报错。
- **API：** `read_fasta_qual_pairs(fasta_path[必须], qual_path[必须], max_records[可选], verbose[可选])`
- **输入：** 本地 FASTA 文件与 QUAL 文件。
- **输出：** 字典，含 `record_count`、`records`（每条含 id、description、sequence、qualities）与 `provenance`。

**示例文件：** `reads.fa`

```text
>r1 desc
ACGT
>r2
GGCC
```

**示例文件：** `reads.qual`

```text
>r1 desc
30 20 10 5
>r2
40 40 40 40
```

**示例代码：**

```python
from dnakit.sequence import read_fasta_qual_pairs

result = read_fasta_qual_pairs("reads.fa", "reads.qual", verbose=False)
print(result["record_count"])
for record in result["records"]:
    print(record["id"], record["sequence"], record["qualities"])
```

**示例输出：**

```text
2
r1 ACGT [30, 20, 10, 5]
r2 GGCC [40, 40, 40, 40]
```

## 10) `write_qual_file`

- **作用：** 把 DNA 记录的 PHRED 质量值写入一个新的 QUAL 文件。
- **处理方法：** 每条记录的质量值数量必须与碱基数相同，取值 0 到 93。使用 Biopython 的 `qual` 格式写出器。
- **API：** `write_qual_file(records[必须], output_path[必须], verbose[可选])`
- **输入：** 记录列表（每条含 `id`、`sequence` 与 `qualities`）以及新输出路径。
- **输出：** 字典，含 `output_path`、`record_count` 与 `provenance`。

**示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import write_qual_file

out = Path(tempfile.mkdtemp()) / "reads.qual"
record = {"id": "r1", "sequence": "ACGT", "qualities": [30, 20, 10, 5]}
result = write_qual_file([record], out, verbose=False)
print(result["record_count"])
print(out.read_text().strip())
```

**示例输出：**

```text
1
>r1
30 20 10 5
```

## 11) `melting_temperature_gc`

- **作用：** 用 GC 含量的经验公式估算无歧义 DNA 引物的熔解温度。
- **处理方法：** 调用 Bio.SeqUtils.MeltingTemp.Tm_GC，参数 check 与 strict 均为真，因此含 N 等歧义碱基会报错。`valueset` 为 1 到 7 的系数组；盐浓度单位为 mM；`saltcorr` 为 0 到 7 的盐校正方法。这不是最近邻模型，最近邻模型见现有 `melting_temperature`。
- **API：** `melting_temperature_gc(sequence[必须], valueset[可选], na[可选], k[可选], tris[可选], mg[可选], dntps[可选], saltcorr[可选], mismatch[可选], verbose[可选])`
- **输入：** 一个无歧义 DNA 字符串与盐、校正参数。
- **输出：** 字典，含 `length`、`tm_celsius`、`valueset`、`saltcorr`、`salt_mM` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import melting_temperature_gc

result = melting_temperature_gc("ACGTACGTACGTACGTAAGC", verbose=False)
print(result["length"], round(result["tm_celsius"], 2))
```

**示例输出：**

```text
20 50.4
```

## 12) `write_gaf_file`

- **作用：** 把基因本体（GO）注释写入一个新的 GAF 2.0 文件。
- **处理方法：** 每条记录是以 GAF 2.0 字段名为键的映射，缺失字段写为空。列表值以 `|` 连接。含制表符、换行或 `|` 的文本会被拒绝。使用 Biopython 的 `writerec`。这是 `writerec` 与 `writebyproteinrec` 的代表函数。
- **API：** `write_gaf_file(records[必须], output_path[必须], verbose[可选])`
- **输入：** GAF 记录列表与新输出路径。
- **输出：** 字典，含 `output_path`、`record_count`、`gaf_version`（固定为 2.0）与 `provenance`。

**示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import write_gaf_file

record = {
    "DB": "UniProtKB",
    "DB_Object_ID": "P12345",
    "DB_Object_Symbol": "GENE1",
    "GO_ID": "GO:0005634",
    "Evidence": "IDA",
    "Aspect": "C",
    "Taxon_ID": ["taxon:9606"],
}
out = Path(tempfile.mkdtemp()) / "annotations.gaf"
result = write_gaf_file([record], out, verbose=False)
fields = out.read_text().rstrip("\n").split("\t")
print(result["record_count"], result["gaf_version"], len(fields))
print(fields[0], fields[1], fields[2], fields[4], fields[12])
```

**示例输出：**

```text
1 2.0 17
UniProtKB P12345 GENE1 GO:0005634 taxon:9606
```

## 13) `read_mast_report`

- **作用：** 读取本地 MEME Suite 的 MAST XML 报告，得到基序、搜索设置与每条序列的命中图示。
- **处理方法：** 调用 Bio.motifs.mast.read。图示按顺序列出命中：数字为间隔长度，`[+1]` 表示基序 1 的正向命中，`[-1]` 表示反向互补命中。格式不正确时报 `ConfigurationError`。
- **API：** `read_mast_report(path[必须], verbose[可选])`
- **输入：** 一个本地 MAST XML 文件。
- **输出：** 字典，含 `version`、`database`、`alphabet`、`strand_handling`、`motifs`、`sequences`（名称与图示）与 `provenance`。

**示例文件：** `mast.xml`

```text
<?xml version="1.0" encoding="UTF-8"?>
<mast version="5.5.0">
<sequence_dbs><sequence_db source="seqs.fa"/></sequence_dbs>
<alphabet name="DNA"/>
<settings strand_handling="combine"/>
<motifs>
<motif id="MOTIF1" alt="m1" length="4"/>
<motif id="MOTIF2" alt="m2" length="3"/>
</motifs>
<sequences>
<sequence name="s1" length="20">
<seg start="1"><hit idx="0" pos="3" rc="n"/></seg>
</sequence>
<sequence name="s2" length="20"/>
</sequences>
</mast>
```

**示例代码：**

```python
from dnakit.sequence import read_mast_report

result = read_mast_report("mast.xml", verbose=False)
print(result["version"], result["alphabet"], result["strand_handling"])
print([motif["name"] for motif in result["motifs"]])
for item in result["sequences"]:
    print(item["name"], item["diagram"])
```

**示例输出：**

```text
5.5.0 DNA combine
['1', '2']
s1 2-[+1]-14
s2 20
```

## 14) `write_motifs`

- **作用：** 把 DNA 基序（每个基序若干等长位点）格式化为 JASPAR、PFM、TRANSFAC 或 ClusterBuster 文本。
- **处理方法：** 对每个基序调用 Bio.motifs.create，再用 Bio.motifs.write 格式化。名称按字母顺序输出；`pfm` 只能包含一个基序。结果返回文本、每个基序的长度与共有序列。
- **API：** `write_motifs(motifs[必须], file_format[必须], verbose[可选])`
- **输入：** 映射：基序名称到位点列表，例如 `{"M1": ["ACGT", "ACGA"]}`，以及格式名。
- **输出：** 字典，含 `file_format`、`motif_count`、`motifs`（名称、长度、共有序列）、`text` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import write_motifs

result = write_motifs({"M1": ["ACGT", "ACGA", "ACGT"]}, "jaspar", verbose=False)
print(result["motif_count"], result["motifs"][0]["consensus"])
print(result["text"].strip())
```

**示例输出：**

```text
1 ACGT
>M1 M1
A [  3.00   0.00   0.00   1.00]
C [  0.00   3.00   0.00   0.00]
G [  0.00   0.00   3.00   0.00]
T [  0.00   0.00   0.00   2.00]
```

## 15) `build_bme_tree`

- **作用：** 用对齐 DNA 的 p-距离构建平衡最小进化树。
- **处理方法：** 先计算 p-距离：不一致位点数除以两条序列都为无歧义碱基的位点数（忽略 gap 与歧义码）。之后调用 skbio.tree.bme。若某对序列没有可比较位点，报错。
- **API：** `build_bme_tree(sequences[必须], neg_as_zero[可选], verbose[可选])`
- **输入：** 映射：分类单元名称到等长的对齐 DNA 字符串（允许 `-`），3 到 500 个分类单元。
- **输出：** 字典，含 `taxa`（按名称排序）、`distance`（固定为 p-distance）、`newick` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import build_bme_tree

sequences = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "TCGTACGTAC", "D": "TCGTTCGTAG"}
result = build_bme_tree(sequences, verbose=False)
print(result["taxa"], result["distance"])
print(result["newick"])
```

**示例输出：**

```text
['A', 'B', 'C', 'D'] p-distance
(B:0.1,(C:0.0,D:0.2):0.1,A:0.0);
```

## 16) `build_gme_tree`

- **作用：** 用对齐 DNA 的 p-距离构建贪心最小进化树。
- **处理方法：** 与 `build_bme_tree` 使用同样的 p-距离，之后调用 skbio.tree.gme。
- **API：** `build_gme_tree(sequences[必须], neg_as_zero[可选], verbose[可选])`
- **输入：** 同 `build_bme_tree`。
- **输出：** 字典，含 `taxa`、`distance`、`newick` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import build_gme_tree

sequences = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "TCGTACGTAC", "D": "TCGTTCGTAG"}
result = build_gme_tree(sequences, verbose=False)
print(result["taxa"], result["distance"])
print(result["newick"])
```

**示例输出：**

```text
['A', 'B', 'C', 'D'] p-distance
(B:0.1,(C:0.0,D:0.2):0.1,A:0.0);
```

## 17) `improve_tree_nni`

- **作用：** 对 Newick 树做最近邻交换（NNI），使其更符合对齐 DNA 的 p-距离。
- **处理方法：** 叶子名称必须与序列名称完全一致，且至少四个叶子（三叶树没有内部边可交换）。调用 skbio.tree.nni；上游的 ValueError 转为 `ConfigurationError`。
- **API：** `improve_tree_nni(newick[必须], sequences[必须], balanced[可选], neg_as_zero[可选], verbose[可选])`
- **输入：** Newick 字符串与对齐 DNA 映射。
- **输出：** 字典，含 `tip_count`、`balanced`、`neg_as_zero`、`distance`、`newick` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import improve_tree_nni

sequences = {"A": "ACGTACGTAC", "B": "ACGTACGAAC", "C": "TCGTACGTAC", "D": "TCGTTCGTAG"}
result = improve_tree_nni("((A:1,C:1):1,(B:1,D:1):1);", sequences, verbose=False)
print(result["tip_count"], result["balanced"])
print(result["newick"])
```

**示例输出：**

```text
4 True
(B:0.1,(C:0.0,D:0.2):0.1,A:0.0);
```

## 18) `tree_path_distances`

- **作用：** 用树中叶子两两之间的路径长度比较多棵树。
- **处理方法：** 调用 skbio.tree.path_dists，精确计算，不做抽样。`metric` 可为 euclidean、cityblock 或 sqeuclidean；`use_length` 决定是否使用分支长度，否则每条边计为 1。
- **API：** `tree_path_distances(newicks[必须], ids[可选], metric[可选], use_length[可选], shared_by_all[可选], verbose[可选])`
- **输入：** 2 到 200 个 Newick 字符串，可选标签列表（唯一、无空白）。
- **输出：** 字典，含 `ids`、`matrix`、`metric`、`use_length`、`shared_by_all` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import tree_path_distances

newicks = ["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"]
result = tree_path_distances(newicks, ids=["t1", "t2"], verbose=False)
print(result["ids"])
print(result["matrix"])
```

**示例输出：**

```text
['t1', 't2']
[[0.0, 4.0], [4.0, 0.0]]
```

## 19) `weighted_robinson_foulds_distances`

- **作用：** 用加权 Robinson-Foulds 距离比较多棵树。
- **处理方法：** 调用 skbio.tree.wrf_dists。每条分支（split）带有分支长度；`metric` 合并共有分支的长度差与独有分支的长度。`include_tips` 决定是否计入末端分支，`rooted` 决定是否按有根树比较。
- **API：** `weighted_robinson_foulds_distances(newicks[必须], ids[可选], metric[可选], rooted[可选], include_tips[可选], shared_by_all[可选], verbose[可选])`
- **输入：** 2 到 200 个 Newick 字符串，可选标签列表。
- **输出：** 字典，含 `ids`、`matrix`、`metric`、`rooted`、`include_tips`、`shared_by_all` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import weighted_robinson_foulds_distances

newicks = ["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"]
result = weighted_robinson_foulds_distances(newicks, ids=["t1", "t2"], verbose=False)
print(result["metric"], result["matrix"])
```

**示例输出：**

```text
cityblock [[0.0, 4.0], [4.0, 0.0]]
```

## 20) `adapter_untrimmed_interval`

- **作用：** 返回一条 reads 在去除一个或多个接头后仍保留的区间。
- **处理方法：** 接头按给定顺序应用，每个接头只在前一个接头之后剩下的部分中匹配，与 cutadapt 的链式接头相同。未保留区间由 `cutadapt.adapters.remainder` 计算；无匹配时整条 reads 保留。
- **API：** `adapter_untrimmed_interval(sequence[必须], adapters[必须], position[可选], max_errors[可选], min_overlap[可选], verbose[可选])`
- **输入：** 一条 DNA 字符串、1 到 10 条接头序列，以及接头位置（`3prime`、`5prime` 或 `anywhere`）。
- **输出：** 字典，含 `start`、`stop`（0 起、半开区间）、`length`、`untrimmed_sequence`、`adapter_matches`（每个匹配的接头、起止与错误数）与 `provenance`。

**示例代码：**

```python
from dnakit.sequencing_qc import adapter_untrimmed_interval

result = adapter_untrimmed_interval("ACGTACGTTTAGATCGGAA", ["AGATCGG"], verbose=False)
print(result["start"], result["stop"], result["untrimmed_sequence"])
match = result["adapter_matches"][0]
print(match["start"], match["stop"], match["errors"])
```

**示例输出：**

```text
0 10 ACGTACGTTT
10 17 0
```

## 21) `barcode_from_bam_tag`

- **作用：** 从单条比对记录的字符串标签中提取 UMI，并可选提取细胞条形码。
- **处理方法：** 在内存中构造 pysam 记录并调用 `umi_tools.sam_methods.get_barcode_tag`。UMI 取自 `umi_tag`（默认 RX），可按 `umi_tag_split` 截断、按 `umi_tag_delim` 删除字符；细胞条形码同理，`cell_tag_split` 默认为 `-`，用于去掉 10x 的 GEM 后缀。缺少所需标签时报错。
- **API：** `barcode_from_bam_tag(tags[必须], umi_tag[可选], umi_tag_split[可选], umi_tag_delim[可选], cell_barcode[可选], cell_tag[可选], cell_tag_split[可选], cell_tag_delim[可选], verbose[可选])`
- **输入：** 映射：两字母 SAM 标签名到字符串值。
- **输出：** 字典，含 `umi`、`cell_barcode`（未请求时为 None）、`umi_tag`、`cell_tag` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequencing_qc import barcode_from_bam_tag

umi = barcode_from_bam_tag({"RX": "ACGTAC-1"}, umi_tag_split="-", verbose=False)
cell = barcode_from_bam_tag(
    {"RX": "ACGTAC-1", "CB": "AAACCT-1"}, cell_barcode=True, cell_tag="CB", verbose=False
)
print(umi["umi"], umi["cell_barcode"])
print(cell["umi"], cell["cell_barcode"])
```

**示例输出：**

```text
ACGTAC None
ACGTAC-1 AAACCT
```

## 22) `metacontig_contig_map`

- **作用：** 把基因（metacontig）映射到本地 SAM 或 BAM 中属于该基因的参考序列。
- **处理方法：** 调用 `umi_tools.sam_methods.getMetaContig2contig`。基因-转录本文件为两列制表符文本；只保留比对头部中存在的转录本。与 umi_tools 一致：`#` 开头的行被跳过，遇到第一个空行即结束。
- **API：** `metacontig_contig_map(bam_path[必须], gene_transcript_path[必须], verbose[可选])`
- **输入：** 本地 SAM 或 BAM 文件，以及基因-转录本两列文件。
- **输出：** 字典，含 `metacontig_count`、`metacontigs`（基因到排序后的参考序列列表）与 `provenance`。

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:coordinate
@SQ	SN:tx1	LN:100
@SQ	SN:tx2	LN:100
@SQ	SN:tx3	LN:100
r1	0	tx1	1	60	5M	*	0	0	ACGTA	IIIII	RX:Z:ACGTAC-1
r2	0	tx2	3	60	4M	*	0	0	ACGT	IIII	RX:Z:GGGGTT-1
r3	16	tx2	10	60	3S5M	*	0	0	ACGTACGT	IIIIIIII
```

**示例文件：** `gene_to_transcript.tsv`

```text
GENE1	tx1
GENE1	tx2
GENE2	tx3
```

**示例代码：**

```python
from dnakit.sequencing_qc import metacontig_contig_map

result = metacontig_contig_map("reads.sam", "gene_to_transcript.tsv", verbose=False)
print(result["metacontig_count"])
for gene, transcripts in result["metacontigs"].items():
    print(gene, transcripts)
```

**示例输出：**

```text
2
GENE1 ['tx1', 'tx2']
GENE2 ['tx3']
```

## 23) `tag_reads_by_metacontig`

- **作用：** 在带索引的 BAM 上按基因（metacontig）顺序统计比对记录数。
- **处理方法：** 调用 `umi_tools.sam_methods.metafetcher`，在内存中为每条记录打上 `metatag` 标签，不写出任何文件。统计按基因进行，共享参考序列的两个基因各计一次。BAM 需要 `.bai` 或 `.csi` 索引。
- **API：** `tag_reads_by_metacontig(bam_path[必须], gene_transcript_path[必须], metatag[可选], verbose[可选])`
- **输入：** 带索引的本地 BAM 文件与基因-转录本文件。
- **输出：** 字典，含 `metacontig_count`、`read_count`、`reads_per_metacontig`（每个基因的记录数）与 `provenance`。

**示例代码：**

```python
import tempfile
from pathlib import Path

import pysam

from dnakit.sequencing_qc import tag_reads_by_metacontig

bam = Path(tempfile.mkdtemp()) / "reads.bam"
with pysam.AlignmentFile("reads.sam", "r") as source:
    with pysam.AlignmentFile(str(bam), "wb", template=source) as target:
        for record in source:
            target.write(record)
pysam.index(str(bam))
result = tag_reads_by_metacontig(bam, "gene_to_transcript.tsv", verbose=False)
print(result["metacontig_count"], result["read_count"], result["reads_per_metacontig"])
```

**示例输出：**

```text
2 3 {'GENE1': 3, 'GENE2': 0}
```

## 24) `soft_clip_adjusted_position`

- **作用：** 根据 CIGAR 计算 UMI 去重所用的位置，并判断是否为剪接读段。
- **处理方法：** 调用 `umi_tools.sam_methods.get_read_position`。正向读段的位置是去掉前端软剪切后的起点；反向读段的位置是加上末端软剪切后的终点。含 `N` 的 CIGAR，或 5' 端软剪切长度超过 `soft_clip_threshold`（默认 4，与 umi_tools 的默认值一致）时，视为剪接。
- **API：** `soft_clip_adjusted_position(cigar[必须], reference_start[可选], is_reverse[可选], soft_clip_threshold[可选], verbose[可选])`
- **输入：** CIGAR 字符串、0 起的比对起点与链方向。
- **输出：** 字典，含 `start`、`position`、`is_spliced`、`splice_offset`（第一个剪接点距读段起点的偏移，无则为 None）与 `provenance`。

**示例代码：**

```python
from dnakit.sequencing_qc import soft_clip_adjusted_position

fwd = soft_clip_adjusted_position("3S5M2N2M", reference_start=10, verbose=False)
rev = soft_clip_adjusted_position("5M2S", reference_start=10, is_reverse=True, verbose=False)
print(fwd["position"], fwd["splice_offset"])
print(rev["start"], rev["position"], rev["is_spliced"])
```

**示例输出：**

```text
7 8
10 17 False
```

## 25) `parse_prefixed_umi_name`

- **作用：** 从 umi_tools extract 写入的读段名称中读取 UMI，并可选读取细胞条形码。
- **处理方法：** 调用 `umi_tools.sam_methods.get_barcode_umis`。名称按 `:` 拆分；`UMI_<序列>` 字段给出 UMI，`CELL_<序列>` 字段给出细胞条形码（需要 `cell_barcode=True`）。没有 UMI 字段时报错。
- **API：** `parse_prefixed_umi_name(read_name[必须], cell_barcode[可选], verbose[可选])`
- **输入：** 一个读段名称字符串。
- **输出：** 字典，含 `umi`、`cell_barcode`、`cell_barcode_requested` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequencing_qc import parse_prefixed_umi_name

result = parse_prefixed_umi_name("HWI:1:UMI_ACGT:CELL_TTAA", cell_barcode=True, verbose=False)
print(result["umi"], result["cell_barcode"])
```

**示例输出：**

```text
ACGT TTAA
```

## 26) `resolve_above_knee_barcodes`

- **作用：** 丢弃（或校正）与更常见白名单条形码仅相差少量编辑的条形码。
- **处理方法：** 调用 `umi_tools.whitelist_methods.errorDetectAboveThreshold`。条形码从计数最低者开始检查，若其与更常见条形码有近似（不超过 `errors` 处编辑）则丢弃。`resolution_method="correct"` 时，单个替换的近似会记入 `true_to_false_map`。输入字典不会被修改。
- **API：** `resolve_above_knee_barcodes(cell_barcode_counts[必须], cell_whitelist[必须], true_to_false_map[必须], errors[可选], resolution_method[可选], verbose[可选])`
- **输入：** 条形码计数映射、白名单列表与真-假条形码映射。
- **输出：** 字典，含 `kept`、`discarded`、`true_to_false_map`、`resolution_method`、`errors` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequencing_qc import resolve_above_knee_barcodes

result = resolve_above_knee_barcodes(
    {"AAAA": 10, "AAAT": 5}, ["AAAA", "AAAT"], {}, verbose=False
)
print(result["kept"], result["discarded"])
```

**示例输出：**

```text
['AAAA'] ['AAAT']
```
