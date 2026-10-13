# H3 族变体复核：Bio、skbio、screed、cutadapt、umi_tools、khmer

本页记录对 108 行 `family-variant` 台账的逐行复核，以及因此新增的 12 个包装函数（10 个来自 Biopython，1 个来自 scikit-bio，1 个来自 khmer）。新函数位于 `dnakit.sequence`、`dnakit.sequence` 和 `dnakit.microbiome`，已从各区域包导出，示例从区域包导入。

- 离线运行：所有函数只读取调用参数或本地文件，不访问网络，也不写文件。
- 判定规则：别名、重导出或同一对象记 `duplicate-of`；只有当现有包装函数的某个参数确实能得到变体行为（用 `inspect.signature` 核对参数存在、实际调用核对结果）才保留 `family-variant`；否则逐一记为 `wrapped`、`deprecated`、`broken-upstream` 或 `internal`。

## 复核结果

| 包 | 原 family-variant | duplicate-of | family-variant（已核实） | wrapped | 其他 |
| --- | ---: | ---: | ---: | ---: | --- |
| Bio | 56 | 9 | 24 | 17 | deprecated 3、broken-upstream 1、internal 2 |
| skbio | 18 | 5 | 12 | 1 | 无 |
| screed | 5 | 3 | 2 | 0 | 无 |
| cutadapt | 13 | 6 | 2 | 0 | internal 5 |
| umi_tools | 8 | 2 | 4 | 0 | internal 2 |
| khmer | 8 | 0 | 5 | 3 | 无 |

主要发现（原声明有误的地方）：

- `Bio.Blast.parse/read` 能读 BLAST XML2，而 `parse_search_results` 背后的 SearchIO 对 XML2 文件不报错却返回空结果；`Bio.Blast.write(fmt='XML2')` 在 Biopython 1.88 中因 `XML2Writer` 缺少 `_write_query_id` 而失败。
- `Bio.motifs.pfm.write` 不是 JASPAR 的 pfm 布局，其输出与 `write_motifs(file_format='clusterbuster')` 相同。
- `Bio.UniProt.GOA.gafbyproteiniterator` 丢掉最后一个蛋白分组（P1,P1,P2,P3 只得到前两组），记 `broken-upstream`；`as_tab`、`as_fastq_illumina`、`as_fastq_solexa` 已被 Biopython 弃用。
- `khmer.load_counttable`：`query_saved_countgraph` 读 counttable 文件会静默返回 0（哈希布局不同），所以不是它的变体；nodetable 与 nodegraph 同理，且 `inspect_khmer_graph` 只读文件头。
- `skbio.sequence.distance.jc69_correct` 输入是 p 距离，`nucleotide_distance(model='jc69')` 无法校正外部 p 距离。
- 多记录读取（多棵树、多个比对块、多条 PHD 或 Primer3 记录、多条 GenBank 头信息）无法通过现有单记录包装函数的参数得到，已新增函数。
- `cutadapt.align.naive/py/slow_edit_environment` 是测试参考实现，`umi_tools.network.recursive_search` 把状态存放在函数属性里，均记 `internal`。

现有 `Bio.Phylo.parse`、`Bio.Blast.parse` 的 `duplicate-of` 声明不在本组行内，但上述发现表明它们并不完全等价（多棵树、XML2）。

## 1) `read_alignment_file`

- **作用：** 读取 BED、Clustal、FASTA、MAF、PHYLIP、PSL 或 SAM 文件中的成对比对。
- **处理方法：** 调用 Bio.Align.parse；`exactly_one=True` 时调用 Bio.Align.read，文件中不是恰好一条比对就抛出 `ConfigurationError`。每条比对给出序列名、比对列数、分值（格式带分值时）、分块坐标（每条序列一行）以及带缺口的行（文件带序列时；PSL、SAM、BED 不带序列，此项为 null）。现有 `read_psl_alignments` 只读 PSL。
- **API：** `read_alignment_file(path[必须], file_format[可选], exactly_one[可选], max_alignments[可选], verbose[可选])`
- **输入：** 文件路径、格式名、是否要求恰好一条比对、比对条数上限。
- **输出：** 字典，含 `alignment_count`、`alignments`（`sequence_ids`、`aligned_columns`、`score`、`coordinates`、`aligned_rows`）与 `provenance`。

**示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.sequence import read_alignment_file
from dnakit.sequence import write_pairwise_alignment_file

path = Path(tempfile.mkdtemp()) / "pair.maf"
write_pairwise_alignment_file("ACGTACGGTT", "ACGTACGTT", path, file_format="maf", verbose=False)
result = read_alignment_file(path, file_format="maf", exactly_one=True, verbose=False)
item = result["alignments"][0]
print(result["alignment_count"], item["sequence_ids"], item["score"])
print(item["aligned_rows"])
print(item["coordinates"])
```

**示例输出：**

```text
1 ['query', 'target'] 16.0
['ACGTACGGTT', 'ACGTAC-GTT']
[[0, 6, 7, 10], [0, 6, 6, 9]]
```

## 2) `read_alignment_blocks`

- **作用：** 读取 MAF 或 FASTA-M10 文件中的全部比对块。
- **处理方法：** 调用 Bio.AlignIO.parse。MAF 文件每个比对区域一块，FASTA-M10 每个命中一块；现有 `read_alignment` 使用 AlignIO.read，遇到多块会报 More than one record。每块给出各行的序列名、带缺口序列和行注释（MAF 的 `start`、`size`、`strand`、`srcSize`），以及块注释。
- **API：** `read_alignment_blocks(path[必须], file_format[可选], max_blocks[可选], verbose[可选])`
- **输入：** 文件路径、格式名（`maf` 或 `fasta-m10`）、块数上限。
- **输出：** 字典，含 `block_count`、`blocks`（`rows`、`alignment_length`、`annotations`）与 `provenance`。

**示例文件：** `blocks.maf`

```text
##maf version=1 scoring=none
a score=23262.0
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG

a score=5062.0
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA

```

**示例代码：**

```python
from dnakit.sequence import read_alignment_blocks

result = read_alignment_blocks("blocks.maf", verbose=False)
print(result["block_count"])
for block in result["blocks"]:
    first = block["rows"][0]
    print(block["row_count"], block["alignment_length"], first["id"], first["annotations"]["start"])
```

**示例输出：**

```text
2
3 42 hg18.chr7 27578828
3 6 hg18.chr7 27699739
```

## 3) `read_substitution_matrix_file`

- **作用：** 读取自定义核苷酸替换矩阵文件。
- **处理方法：** 调用 Bio.Align.substitution_matrices.read。文件含可选的 `#` 注释行、符号表头行和每个符号一行。字母表只能含 A、C、G、T、U 与 IUPAC 符号，否则抛出 `ConfigurationError`。内置矩阵仍用 `load_substitution_matrix` 按名字加载。
- **API：** `read_substitution_matrix_file(path[必须], verbose[可选])`
- **输入：** 矩阵文件路径。
- **输出：** 字典，含 `alphabet`、`size`、`symmetric`、`matrix`、`comments` 与 `provenance`。

**示例文件：** `nuc.mat`

```text
# demo nucleotide matrix
   A  C  G  T
A  2 -3 -3 -3
C -3  2 -3 -3
G -3 -3  2 -3
T -3 -3 -3  2
```

**示例代码：**

```python
from dnakit.sequence import read_substitution_matrix_file

result = read_substitution_matrix_file("nuc.mat", verbose=False)
print(result["alphabet"], result["symmetric"], result["comments"])
print(result["matrix"][0])
```

**示例输出：**

```text
ACGT True ['demo nucleotide matrix']
[2.0, -3.0, -3.0, -3.0]
```

## 4) `read_blast_records`

- **作用：** 读取 BLAST XML 或 XML2 输出，每个查询一条记录。
- **处理方法：** 调用 Bio.Blast.parse；`exactly_one=True` 时调用 Bio.Blast.read。每条记录给出查询、程序、版本、数据库，以及每个命中的高分片段对（分值、bit score、E 值、一致数、靶标与查询的坐标，0 起始半开）。`parse_search_results` 背后的 SearchIO 解析器读取 XML2 布局时不报错但返回空结果。
- **API：** `read_blast_records(path[必须], exactly_one[可选], max_records[可选], max_hits[可选], verbose[可选])`
- **输入：** 文件路径、是否要求恰好一个查询、记录数与命中数上限。
- **输出：** 字典，含 `record_count`、`records`（`hits`、`hsps`）与 `provenance`。

**示例文件：** `blast_small.xml`

```text
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.13.0+</BlastOutput_version>
  <BlastOutput_db>demo_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>query_one demo query</BlastOutput_query-def>
  <BlastOutput_query-len>20</BlastOutput_query-len>
    <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>query_one demo query</Iteration_query-def>
      <Iteration_query-len>20</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>subject_one demo subject</Hit_def>
          <Hit_accession>0</Hit_accession>
          <Hit_len>60</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>37.4</Hsp_bit-score>
              <Hsp_score>20</Hsp_score>
              <Hsp_evalue>1.2e-05</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>11</Hsp_hit-from>
              <Hsp_hit-to>30</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>20</Hsp_identity>
              <Hsp_positive>20</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>20</Hsp_align-len>
              <Hsp_qseq>ACGTACGGTTAACCGGATAC</Hsp_qseq>
              <Hsp_hseq>ACGTACGGTTAACCGGATAC</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
          </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
```

**示例代码：**

```python
from dnakit.sequence import read_blast_records

result = read_blast_records("blast_small.xml", verbose=False)
record = result["records"][0]
hit = record["hits"][0]
hsp = hit["hsps"][0]
print(record["query_id"], record["query_length"], record["program"], record["database"])
print(hit["target_id"], hit["target_length"])
print(hsp["evalue"], hsp["bit_score"], hsp["target_coordinates"], hsp["query_coordinates"])
```

**示例输出：**

```text
Query_1 20 blastn demo_db
gnl|BL_ORD_ID|0 60
1.2e-05 37.4 [10, 30] [0, 20]
```

## 5) `parse_primer3_records`

- **作用：** 读取 EMBOSS eprimer3 或 Primer3 报告中的每一条记录（每个靶序列一条）。
- **处理方法：** 调用 Bio.Emboss.Primer3.parse。现有 `parse_primer3_output` 用 Primer3.read，遇到多靶序列报告会抛出 `ConfigurationError`。每条记录给出头部注释以及引物对（正向、反向、可选内部寡核苷酸）。空文件或没有记录的文件抛出 `ConfigurationError`（上游 parse 对空流会死循环，这里提前拦截）。
- **API：** `parse_primer3_records(path[必须], max_records[可选], verbose[可选])`
- **输入：** 文件路径与记录数上限。
- **输出：** 字典，含 `record_count`、`pair_count`、`records` 与 `provenance`。

**示例文件：** `primer3_multi.txt`

```text
# EPRIMER3 RELEASE 5.0.0
#                      Start  Len   Tm     GC%   Sequence

   1 PRODUCT SIZE: 60
     FORWARD PRIMER    1  20  60.00  50.00  ACGTACGTACGTACGTACGT
     REVERSE PRIMER   60  20  60.00  50.00  TGCATGCATGCATGCATGCA

# EPRIMER3 RELEASE 5.0.0
#                      Start  Len   Tm     GC%   Sequence

   1 PRODUCT SIZE: 80
     FORWARD PRIMER    5  18  59.50  55.00  GGCCTTAAGGCCTTAAGG
     REVERSE PRIMER   80  18  59.80  50.00  CCGGAATTCCGGAATTCC
```

**示例代码：**

```python
from dnakit.sequence import parse_primer3_records

result = parse_primer3_records("primer3_multi.txt", verbose=False)
print(result["record_count"], result["pair_count"])
for record in result["records"]:
    pair = record["primer_pairs"][0]
    print(pair["product_size"], pair["forward"]["sequence"], pair["reverse"]["start"])
```

**示例输出：**

```text
2 2
60 ACGTACGTACGTACGTACGT 60
80 GGCCTTAAGGCCTTAAGG 80
```

## 6) `read_genbank_record_fields`

- **作用：** 读取 GenBank 记录的头部字段和特征表。
- **处理方法：** 调用 Bio.GenBank.parse；`exactly_one=True` 时调用 Bio.GenBank.read。底层解析器保留 `read_sequence_records` 丢弃的字段：LOCUS 行（名称、长度、分子类型、分类部门、日期）、登录号、版本、关键词、来源与物种、分类谱系、参考文献数、注释，以及每个特征的位置文本和限定词。
- **API：** `read_genbank_record_fields(path[必须], exactly_one[可选], max_records[可选], max_bases[可选], verbose[可选])`
- **输入：** 文件路径、是否要求恰好一条记录、记录数与碱基数上限。
- **输出：** 字典，含 `record_count`、`total_bases`、`records` 与 `provenance`。

**示例文件：** `two.gb`

```text
LOCUS       SEQ1                      18 bp    DNA     linear   BCT 01-JAN-2020
DEFINITION  Demo sequence one.
ACCESSION   SEQ1
VERSION     SEQ1.1
KEYWORDS    demo; test.
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria; Proteobacteria; Gammaproteobacteria.
FEATURES             Location/Qualifiers
     source          1..18
                     /organism="Escherichia coli"
     gene            1..12
                     /gene="abcA"
ORIGIN
        1 atgaaacgta ttgcgtaa
//
LOCUS       SEQ2                      12 bp    DNA     linear   BCT 02-JAN-2020
DEFINITION  Demo sequence two.
ACCESSION   SEQ2
VERSION     SEQ2.1
KEYWORDS    .
SOURCE      Bacillus subtilis
  ORGANISM  Bacillus subtilis
            Bacteria; Firmicutes; Bacilli.
FEATURES             Location/Qualifiers
     CDS             1..12
                     /gene="xyzB"
ORIGIN
        1 atggcgtacg tt
//
```

**示例代码：**

```python
from dnakit.sequence import read_genbank_record_fields

result = read_genbank_record_fields("two.gb", verbose=False)
print(result["record_count"], result["total_bases"])
for record in result["records"]:
    print(record["locus"], record["data_file_division"], record["date"], record["organism"])
    print(record["taxonomy"], [f["type"] for f in record["features"]])
```

**示例输出：**

```text
2 30
SEQ1 BCT 01-JAN-2020 Escherichia coli
['Bacteria', 'Proteobacteria', 'Gammaproteobacteria'] ['source', 'gene']
SEQ2 BCT 02-JAN-2020 Bacillus subtilis
['Bacteria', 'Firmicutes', 'Bacilli'] ['CDS']
```

## 7) `read_phylogeny_collection`

- **作用：** 读取 Newick、NEXUS、NeXML、phyloXML 或 CDAO 文件中的全部树。
- **处理方法：** 调用 Bio.Phylo.parse，它进入 NexusIO、NeXMLIO、PhyloXMLIO、CDAOIO 的解析器。`read_tree` 只读一棵树，遇到多棵树抛出错误，所以 bootstrap 或后验树集合需要本函数。每棵树给出名称、是否有根、排序后的叶名、总枝长和 Newick 字符串。没有树时抛出 `ConfigurationError`。
- **API：** `read_phylogeny_collection(path[必须], file_format[可选], max_trees[可选], verbose[可选])`
- **输入：** 文件路径、格式名、树数上限。
- **输出：** 字典，含 `tree_count`、`trees` 与 `provenance`。

**示例文件：** `trees.nexus`

```text
#NEXUS
BEGIN TAXA;
DIMENSIONS NTAX=4;
TAXLABELS A B C D;
END;
BEGIN TREES;
TREE t1 = ((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06);
TREE t2 = ((A:0.1,C:0.2):0.05,(B:0.3,D:0.4):0.06);
END;
```

**示例代码：**

```python
from dnakit.sequence import read_phylogeny_collection

result = read_phylogeny_collection("trees.nexus", file_format="nexus", verbose=False)
print(result["tree_count"])
for tree in result["trees"]:
    print(tree["name"], tree["tips"], round(tree["total_branch_length"], 2))
    print(tree["newick"])
```

**示例输出：**

```text
2
t1 ['A', 'B', 'C', 'D'] 1.11
((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.06):0;
t2 ['A', 'B', 'C', 'D'] 1.11
((A:0.1,C:0.2):0.05,(B:0.3,D:0.4):0.06):0;
```

## 8) `bootstrap_replicate_trees`

- **作用：** 生成比对序列集的 bootstrap 重复树。
- **处理方法：** 调用 Bio.Phylo.Consensus.bootstrap_trees：对比对列有放回重抽样，每次重复构建一棵树（`nj` 或 `upgma`）。`bootstrap_consensus_tree` 只返回多数规则一致树；本函数返回重复树本身，可交给 `bootstrap_branch_support` 计数。调用期间用 `seed` 设置全局随机状态，结束后恢复。
- **API：** `bootstrap_replicate_trees(sequences[必须], replicates[可选], method[可选], distance_model[可选], seed[可选], max_length[可选], verbose[可选])`
- **输入：** 等长的 DNA 比对（3 至 200 条）、重复次数、方法、距离模型与种子。
- **输出：** 字典，含 `trees`（Newick 列表）、`ids`、`replicates`、`seed` 与 `provenance`。

**示例代码：**

```python
import random

from dnakit.sequence import bootstrap_replicate_trees

alignment = {
    "a": "ACGTACGTAC",
    "b": "ACGTACGAAC",
    "c": "ACGAACGTTC",
    "d": "TCGTACGTAC",
}
random.seed(5)
expected = random.random()
random.seed(5)
first = bootstrap_replicate_trees(alignment, replicates=4, seed=7, verbose=False)
second = bootstrap_replicate_trees(alignment, replicates=4, seed=7, verbose=False)
print(len(first["trees"]), first["trees"] == second["trees"], first["ids"])
print(all(name in tree for tree in first["trees"] for name in "abcd"))
print(random.random() == expected)
```

**示例输出：**

```text
4 True ['a', 'b', 'c', 'd']
True
True
```

## 9) `read_phd_records`

- **作用：** 读取 PHD 文件中的每一条 PHRED 碱基调用记录。
- **处理方法：** 调用 Bio.Sequencing.Phd.parse。`read_phd_sequence` 用 Phd.read，只接受恰好一条记录。每条记录给出文件名、碱基序列、每个碱基的 PHRED 质量、峰位置（文件没有时为 null）和注释字段。
- **API：** `read_phd_records(path[必须], max_records[可选], verbose[可选])`
- **输入：** 文件路径与记录数上限。
- **输出：** 字典，含 `record_count`、`records` 与 `provenance`。

**示例文件：** `reads.phd`

```text
BEGIN_SEQUENCE read1
BEGIN_COMMENT
CHROMAT_FILE: read1.ab1
PHRED_VERSION: 0.980528.c
END_COMMENT
BEGIN_DNA
A 10 1
C 20 4
END_DNA
END_SEQUENCE

BEGIN_SEQUENCE read2
BEGIN_COMMENT
CHROMAT_FILE: read2.ab1
PHRED_VERSION: 0.980528.c
END_COMMENT
BEGIN_DNA
G 30 2
T 40 3
A 50 5
END_DNA
END_SEQUENCE
```

**示例代码：**

```python
from dnakit.sequence import read_phd_records

result = read_phd_records("reads.phd", verbose=False)
print(result["record_count"])
for record in result["records"]:
    print(record["file_name"], record["sequence"], record["qualities"], record["peak_locations"])
```

**示例输出：**

```text
2
read1 AC [10, 20] [1, 4]
read2 GTA [30, 40, 50] [2, 3, 5]
```

## 10) `read_goa_gene_product_file`

- **作用：** 读取 GOA 基因产物关联（GPA）或基因产物信息（GPI）文件。
- **处理方法：** `file_format='gpa'` 调用 Bio.UniProt.GOA.gpa_iterator（1.0、1.1 版：每行一个基因产物、GO 术语、证据与参考文献）；`'gpi'` 调用 gpi_iterator（1.0、1.1、1.2 版：符号、名称、同义名、类型、分类单元与交叉引用）。版本来自必须存在的 `!gpa-version` 或 `!gpi-version` 头行。`read_gaf_annotations` 只读 GAF 布局。
- **API：** `read_goa_gene_product_file(path[必须], file_format[可选], max_records[可选], verbose[可选])`
- **输入：** 文件路径、格式名（`gpa` 或 `gpi`）、记录数上限。
- **输出：** 字典，含 `version`、`record_count`、`records`（字段名为 GOA 规范名称，列表字段为列表）与 `provenance`。

**示例文件：** `prod.gpa`

```text
!gpa-version: 1.1
UniProtKB	P12345	enables	GO:0003700	PMID:1234567	ECO:0000314			20240101	UniProt		
UniProtKB	P67890	involved_in	GO:0006355	PMID:7654321|PMID:1	ECO:0000501	GO:1		20240102	UniProt		
```

**示例代码：**

```python
from dnakit.sequence import read_goa_gene_product_file

result = read_goa_gene_product_file("prod.gpa", verbose=False)
print(result["version"], result["record_count"])
for record in result["records"]:
    print(record["DB_Object_ID"], record["GO_ID"], record["Qualifier"], record["DB:Reference"], record["With"])
```

**示例输出：**

```text
1.1 2
P12345 GO:0003700 ['enables'] ['PMID:1234567'] []
P67890 GO:0006355 ['involved_in'] ['PMID:7654321', 'PMID:1'] ['GO:1']
```

## 11) `correct_jc69_distance`

- **作用：** 对观测到的 p 距离做 Jukes-Cantor（JC69）校正。
- **处理方法：** 调用 skbio.sequence.distance.jc69_correct：D = -a ln(1 - p/a)，a = (n-1)/n，n 为字母表大小（DNA 为 4）。每个 p 必须在 [0, 1]。p 不小于 a 时（饱和）没有有限距离，`corrected` 中为 null，状态为 `partial`。`nucleotide_distance(model='jc69')` 只能从两条序列出发，无法校正外部测得的 p 距离。
- **API：** `correct_jc69_distance(p_distances[必须], alphabet_size[可选], verbose[可选])`
- **输入：** p 距离列表与字母表大小。
- **输出：** 字典，含 `corrected`、`saturated_count`、`saturation_limit`、`status` 与 `provenance`。

**示例代码：**

```python
from dnakit.sequence import correct_jc69_distance

result = correct_jc69_distance([0.0, 0.1, 0.5, 0.8], verbose=False)
print(result["status"], result["saturation_limit"], result["saturated_count"])
print([None if v is None else round(v, 6) for v in result["corrected"]])
```

**示例输出：**

```text
partial 0.75 1
[0.0, 0.107326, 0.823959, None]
```

## 12) `query_saved_khmer_table`

- **作用：** 查询已保存的 khmer counttable、nodetable 或 nodegraph 中的 k-mer。
- **处理方法：** `table_type` 必须与写出文件的类一致，因为 khmer 不在文件中记录它：`counttable` 调用 khmer.load_counttable 返回计数；`nodetable`、`nodegraph` 分别调用 load_nodetable、load_nodegraph，k-mer 已加入返回 1，否则 0。`query_saved_countgraph` 把 counttable 文件当 countgraph 读时哈希布局不同，会静默返回 0。k-mer 与其反向互补共享一个条目；结果按 k-mer 排序。
- **API：** `query_saved_khmer_table(source[必须], kmers[必须], table_type[可选], verbose[可选])`
- **输入：** 文件路径、k-mer 列表（长度须等于表的 ksize）、表类型。
- **输出：** 字典，含 `values`、`value_kind`（`count` 或 `presence`）、`ksize`、`hashsizes` 与 `provenance`。

**示例代码：**

```python
import tempfile
from pathlib import Path

import khmer

from dnakit.microbiome import query_saved_khmer_table
from dnakit.microbiome import query_saved_countgraph

table = khmer.Counttable(5, 1000, 4)
table.consume("ACGTACGGATCGATCGGCTAGCTAGG")
table.consume("ACGTA")
path = Path(tempfile.mkdtemp()) / "demo.ct"
table.save(str(path))
result = query_saved_khmer_table(path, ["ACGTA", "TACGT", "GGGGG"], verbose=False)
print(result["value_kind"], result["ksize"], result["values"])
print(query_saved_countgraph(path, ["ACGTA"], verbose=False)["counts"])
```

**示例输出：**

```text
count 5 {'ACGTA': 2, 'GGGGG': 0, 'TACGT': 2}
{'ACGTA': 0}
```
