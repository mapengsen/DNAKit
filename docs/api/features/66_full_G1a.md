# 可选序列函数：Biopython 与 scikit-bio 的 DNA 与基因级封装

本页介绍 `dnakit.sequence` 与 `dnakit.sequence` 中的 26 个函数：Biopython 的 20 个，scikit-bio 的 6 个。它们调用上游包中与 DNA 或基因相关的功能，只处理显式给出的字符串、映射或本地文件；不访问网络，也不调用外部程序。当前版本需从这两个模块直接导入，例如 `from dnakit.sequence import read_gaf_annotations`。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。其余上游函数（例如蛋白质结构、生态多样性、远程检索，以及与现有 dnakit 函数重复的实现）未封装，其去向记录在 ledger 中。
## 示例数据

下面的示例使用这些小文件。它们只用于演示接口，数值不代表真实数据。其余示例的输入直接写在代码中。

**示例文件：** `primersearch.txt`

```text
Primer name P1
Amplimer 1
	Sequence: ATGCGATTTGCAAGCTTAGC
	Amplimer length: 42 bp
	Forward primer: 0 mismatches
	Reverse primer: 0 mismatches
```

**示例文件：** `contig.ace`

```text
AS 1 2

CO Contig1 12 2 4 U
ACGTACGTACGT

BQ
30 30 30 30 30 30 30 30 30 30 30 30

AF read1 U 1
AF read2 U 2

RD read1 12 0 0
ACGTACGTACGT

QA 1 12 1 12
DS CHROMAT_FILE read1 PHD_FILE read1.phd.1 TIME Thu Jan 1 00:00:00 2024

RD read2 12 0 0
ACGTACGTACGT

QA 1 12 1 12
DS CHROMAT_FILE read2 PHD_FILE read2.phd.1 TIME Thu Jan 1 00:00:00 2024

```

**示例文件：** `read.phd`

```text
BEGIN_SEQUENCE read1
BEGIN_COMMENT
CHROMAT_FILE: read1.ab1
PHRED_VERSION: 0.980528.c
CALL_METHOD: phred
END_COMMENT
BEGIN_DNA
A 10 1
C 20 4
G 30 2
T 40 3
END_DNA
END_SEQUENCE
```

**示例文件：** `pop.genepop`

```text
Example GenePop file
Locus1
Locus2
Pop
ind1 , 0101 0203
ind2 , 0102 0000
Pop
ind3 , 0201 0303
ind4 , 0000 0303
```

**示例文件：** `annotations.gaf`

```text
!gaf-version: 2.2
UniProtKB	P12345	GENEA		GO:0003700	PMID:1234567	IDA		F	TF A	GENEA|TFA	protein	taxon:9606	20240101	UniProt
UniProtKB	P67890	GENEB	involved_in	GO:0006355	PMID:7654321	IEA		P	Protein B	GENEB	protein	taxon:9606	20240102	UniProt
```

**示例文件：** `gene.txt`

```text
ENTRY       hsa:7157                    CDS       T01001
SYMBOL      TP53
NAME        (RefSeq) tumor protein p53
ORTHOLOGY   K04442  tumor protein p53
ORGANISM    hsa  Homo sapiens (human)
POSITION    17:7661779..7687538
DBLINKS     NCBI-GeneID: 7157
///
```

**示例文件：** `alignments.psl`

```text
psLayout version 3

match	mis-	rep.	N's	Q gap	Q gap	T gap	T gap	strand	Q name	Q size	Q start	Q end	T name	T size	T start	T end	block count	blockSizes	qStarts	tStarts
     	match	match	   	count	bases	count	bases	      									
---------------------------------------------------------------------------------------------------------------
20	0	0	0	0	0	0	0	+	read1	20	0	20	chr1	1000	100	120	1	20,	0,	100,
```

**示例文件：** `primer3.txt`

```text
# PRIMER3 test
#                      Start  Len   Tm     GC%   Sequence
     PRODUCT SIZE: 60
     FORWARD PRIMER    1  20  60.00  50.00  ACGTACGTACGTACGTACGT
     REVERSE PRIMER   60  20  60.00  50.00  TGCATGCATGCATGCATGCA
```

**示例文件：** `reads.fq`

```text
@r1
ACGT
+
IIII
```

## 1) `dnakit.sequence.gc_content_by_codon_position()`

- **作用：** 计算整条序列的 G+C 含量，并分别计算密码子第 1、2、3 位的 G+C 百分比。
- **处理方法：** 调用 Biopython 的 `Bio.SeqUtils.GC123`。输入先去除空白并转为大写；IUPAC 模糊碱基计入长度，但不计为 G 或 C。
- **API：** `dnakit.sequence.gc_content_by_codon_position(sequence[required], max_length[optional], verbose[optional])`。
- **输入：** 一条 DNA 字符串，长度不超过 `max_length`（默认 1,000,000）。
- **输出：** `gc_percent`（包含 `total`、`first`、`second`、`third`，单位为百分比）、`sequence_length` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import gc_content_by_codon_position

result = gc_content_by_codon_position("ATGCGTACGTTAGC", verbose=False)
print(result["gc_percent"])
```

- **示例输出：**

```text
{'total': 50.0, 'first': 40.0, 'second': 60.0, 'third': 50.0}
```

## 2) `dnakit.sequence.gc_skew_windows()`

- **作用：** 按连续、不重叠的窗口计算 GC 偏斜 (G-C)/(G+C)，用于查看复制起点或链偏向。
- **处理方法：** 调用 `Bio.SeqUtils.GC_skew`。最后一个窗口可能较短；不含 G 或 C 的窗口偏斜记为 0.0。窗口起止位置为 1-based 闭区间。
- **API：** `dnakit.sequence.gc_skew_windows(sequence[required], window[optional], max_length[optional], verbose[optional])`。
- **输入：** 一条 DNA 字符串；`window` 为正整数，且不超过序列长度。
- **输出：** `window`、`window_count` 和 `windows`；每个窗口包含 `start`、`end` 与 `skew`。
- **示例代码：**

```python
from dnakit.sequence import gc_skew_windows

result = gc_skew_windows("GGGGCCCCAAAATTTTGGGG", window=5, verbose=False)
for item in result["windows"]:
    print(item["start"], item["end"], item["skew"])
```

- **示例输出：**

```text
1 5 0.6
6 10 -1.0
11 15 0.0
16 20 1.0
```

## 3) `dnakit.sequence.local_composition_complexity()`

- **作用：** 计算每个滑动窗口的局部组成复杂度（LCC），值越低表示组成越偏向少数碱基。
- **处理方法：** 调用 `Bio.SeqUtils.lcc.lcc_mult`，窗口长度由 `window` 指定。输入只允许 A、C、G、T。`values[i]` 对应从第 i+1 位开始的窗口。
- **API：** `dnakit.sequence.local_composition_complexity(sequence[required], window[optional], max_length[optional], verbose[optional])`。
- **输入：** 一条 ACGT 字符串；`window` 不超过序列长度。
- **输出：** `value_count`、`values`、`minimum`、`maximum`、`mean` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import local_composition_complexity

result = local_composition_complexity("ACGTACGTAAAAAAAACGT", window=5, verbose=False)
print(result["value_count"], [round(value, 4) for value in result["values"][:4]])
```

- **示例输出：**

```text
15 [0.961, 0.961, 0.961, 0.961]
```

## 4) `dnakit.sequence.ambiguous_codon_amino_acids()`

- **作用：** 列出一个含 IUPAC 模糊碱基的密码子可能编码的全部氨基酸；只能编码终止子时单独标记。
- **处理方法：** 调用 `Bio.Data.CodonTable.list_possible_proteins`，并使用 NCBI 遗传密码表（默认表 1）。混合了氨基酸和终止子的模糊密码子会报错。
- **API：** `dnakit.sequence.ambiguous_codon_amino_acids(codon[required], genetic_code[optional], verbose[optional])`。
- **输入：** 三个 IUPAC DNA 符号组成的字符串；`genetic_code` 为 1 到 33 之间的整数。
- **输出：** `amino_acids`（排序列表）、`amino_acid_count`、`stop_only` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import ambiguous_codon_amino_acids

for codon in ["GCN", "TTY", "TAR"]:
    result = ambiguous_codon_amino_acids(codon, verbose=False)
    print(codon, result["amino_acids"], result["stop_only"])
```

- **示例输出：**

```text
GCN ['A'] False
TTY ['F'] False
TAR [] True
```

## 5) `dnakit.sequence.expand_ambiguous_codons()`

- **作用：** 从一组具体密码子推出可以无歧义覆盖的 IUPAC 模糊密码子，例如由 TAA 与 TAG 得到 TAR。
- **处理方法：** 调用 `Bio.Data.CodonTable.list_ambiguous_codons`。仅当某个模糊密码子所代表的全部具体密码子都已存在时才会加入，因此 TRR 等不会被加入。
- **API：** `dnakit.sequence.expand_ambiguous_codons(codons[required], max_codons[optional], verbose[optional])`。
- **输入：** 三碱基 ACGT 密码子的序列；重复项会被去除。
- **输出：** `input_codons`、`added_codons`、`expanded_codons`（均已排序）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import expand_ambiguous_codons

result = expand_ambiguous_codons(["TGA", "TAA", "TAG"], verbose=False)
print(result["added_codons"])
print(result["expanded_codons"])
```

- **示例输出：**

```text
['TAR', 'TRA']
['TAA', 'TAG', 'TAR', 'TGA', 'TRA']
```

## 6) `dnakit.sequence.sequence_checksum()`

- **作用：** 计算序列的 SEGUID、CRC32、CRC64 或 GCG 校验值，用于记录和比对序列标识。
- **处理方法：** 调用 `Bio.SeqUtils.CheckSum` 中对应的函数。计算前去除空白并转为大写。SEGUID 与 CRC64 返回字符串，CRC32 与 GCG 返回整数。
- **API：** `dnakit.sequence.sequence_checksum(sequence[required], method[optional], max_length[optional], verbose[optional])`。
- **输入：** 一条 IUPAC DNA 字符串；`method` 为 `seguid`、`crc32`、`crc64` 或 `gcg`。
- **输出：** `method`、`sequence_length`、`value` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import sequence_checksum

for method in ["seguid", "crc32", "crc64", "gcg"]:
    print(method, sequence_checksum("ATGCGT", method=method, verbose=False)["value"])
```

- **示例输出：**

```text
seguid Yrhkea/OcXcnrtHsYwY41PbWdEk
crc32 527699539
crc64 CRC-71A87EB861ADB000
gcg 1573
```

## 7) `dnakit.sequence.mcdonald_kreitman_test()`

- **作用：** 用 McDonald-Kreitman 检验比较两个物种内的多态性与种间固定差异，判断非同义替换是否偏离中性。
- **处理方法：** 调用 `Bio.Align.analysis.mktest`，按密码子计数同义和非同义变化，并以 G 检验给出 p 值。每个 2x2 表格单元都必须非零，否则返回配置错误。
- **API：** `dnakit.sequence.mcdonald_kreitman_test(sequences[required], species[required], genetic_code[optional], max_length[optional], verbose[optional])`。
- **输入：** `sequences` 为 ID 到编码序列的映射，`species` 为 ID 到物种标签的映射。序列需等长且长度为 3 的倍数，只允许 ACGT；至少两个物种。
- **输出：** `p_value`、`sequence_count`、`codon_count`、`species`、`genetic_code` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import mcdonald_kreitman_test

sequences = {"a1": "GCTGAAAAATTT", "a2": "GCCGAAAAATTA",
             "b1": "GCTGAGCAATTT", "b2": "GCTGAGCAATTT"}
species = {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
result = mcdonald_kreitman_test(sequences, species, verbose=False)
print(result["species"], result["codon_count"], round(result["p_value"], 4))
```

- **示例输出：**

```text
['A', 'B'] 4 1.0
```

## 8) `dnakit.sequence.pairwise_dn_ds_matrix()`

- **作用：** 对一组编码序列计算两两 dN（非同义替换率）与 dS（同义替换率）矩阵。
- **处理方法：** 调用 `Bio.Align.analysis.calculate_dn_ds_matrix`，方法为 NG86 或 LWL85。Biopython 用 -1 表示饱和，本函数将其转为 `null` 并计入 `saturated_pair_count`。
- **API：** `dnakit.sequence.pairwise_dn_ds_matrix(sequences[required], method[optional], genetic_code[optional], max_length[optional], verbose[optional])`。
- **输入：** ID 到编码序列的映射；序列需等长，长度为 3 的倍数，只允许 ACGT，至少两条。
- **输出：** `ids`、`dn`、`ds`（对称的完整矩阵，对角线为 0）、`saturated_pair_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import pairwise_dn_ds_matrix

sequences = {
    "s1": "ATGGCTGCCAAATTTGGGCTGAAAGATCGTAAGCCTTTCGCA",
    "s2": "ATGGCAGCCAAATTCGGGCTGAAAGATCGCAAGCCTTTCGCA",
    "s3": "ATGGCTGCCAAATTTGGGCTGAAAGATCGTAAGCCTTTCGCT",
}
result = pairwise_dn_ds_matrix(sequences, verbose=False)
print(result["ids"])
print([round(value, 4) for value in result["ds"][0]])
print(result["saturated_pair_count"])
```

- **示例输出：**

```text
['s1', 's2', 's3']
[0.0, 0.4197, 0.1156]
0
```

## 9) `dnakit.sequence.build_codon_alignment()`

- **作用：** 把 DNA 密码子按蛋白质比对的列结构排列，得到密码子比对。
- **处理方法：** 调用 `Bio.codonalign.build`。蛋白质中的空位会变成三碱基的密码子空位；DNA 的每个密码子必须与对应氨基酸一致。
- **API：** `dnakit.sequence.build_codon_alignment(protein_alignment[required], nucleotide_sequences[required], genetic_code[optional], max_length[optional], verbose[optional])`。
- **输入：** `protein_alignment` 为等长的蛋白质比对（允许 `-`），`nucleotide_sequences` 为无空位的 ACGT 序列；两者的 ID 必须一致。
- **输出：** `ids`、`alignment_length`、`sequences`（ID 到密码子比对字符串）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import build_codon_alignment

protein = {"s1": "MAA-KF", "s2": "MAGKF-", "s3": "MAAKY-"}
nucleotides = {"s1": "ATGGCTGCCAAATTT", "s2": "ATGGCAGGCAAATTT", "s3": "ATGGCTGCCAAATAT"}
result = build_codon_alignment(protein, nucleotides, verbose=False)
print(result["alignment_length"])
for identifier, row in result["sequences"].items():
    print(identifier, row)
```

- **示例输出：**

```text
18
s1 ATGGCTGCC---AAATTT
s2 ATGGCAGGCAAATTT---
s3 ATGGCTGCCAAATAT---
```

## 10) `dnakit.sequence.parse_primersearch_output()`

- **作用：** 读取 EMBOSS primersearch 的结果文件，得到每个引物对应的扩增子序列、长度和附加信息。
- **处理方法：** 调用 `Bio.Emboss.PrimerSearch.read`。本函数只解析已保存的报告，不运行 EMBOSS。引物名称按字典序排列。
- **API：** `dnakit.sequence.parse_primersearch_output(path[required], verbose[optional])`。
- **输入：** primersearch 输出文件的路径。
- **输出：** `primer_count`、`amplicon_count`、`primers`（引物名到扩增子列表）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import parse_primersearch_output

result = parse_primersearch_output("primersearch.txt", verbose=False)
amplicon = result["primers"]["P1"][0]
print(amplicon["sequence"], amplicon["length"])
print(amplicon["details"])
```

- **示例输出：**

```text
ATGCGATTTGCAAGCTTAGC 42
['Forward primer: 0 mismatches', 'Reverse primer: 0 mismatches']
```

## 11) `dnakit.sequence.read_ace_assembly()`

- **作用：** 读取 ACE 组装文件，得到 contig、contig 中的 reads 以及 reads 在 contig 上的位置。
- **处理方法：** 调用 `Bio.Sequencing.Ace.read`。读段的质量剪切与比对剪切位置保持文件中的 1-based 坐标。
- **API：** `dnakit.sequence.read_ace_assembly(path[required], verbose[optional])`。
- **输入：** ACE 文件的路径。
- **输出：** `contig_count`、`read_count`、`contigs`（每个 contig 含 `name`、`length`、`sequence`、`reads` 与 `placements`）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import read_ace_assembly

result = read_ace_assembly("contig.ace", verbose=False)
contig = result["contigs"][0]
print(contig["name"], contig["length"], result["read_count"])
print(contig["placements"])
```

- **示例输出：**

```text
Contig1 12 2
[{'read': 'read1', 'strand': 'U', 'padded_start': 1}, {'read': 'read2', 'strand': 'U', 'padded_start': 2}]
```

## 12) `dnakit.sequence.read_phd_sequence()`

- **作用：** 读取一条 PHRED 碱基识别（PHD）记录，得到碱基序列、逐碱基质量值和峰位置。
- **处理方法：** 调用 `Bio.Sequencing.Phd.read`。没有峰位置的行记为 `null`；注释中的空值会被删除。
- **API：** `dnakit.sequence.read_phd_sequence(path[required], verbose[optional])`。
- **输入：** 一个只包含一条记录的 PHD 文件路径。
- **输出：** `file_name`、`sequence`、`length`、`qualities`、`peak_locations`、`comments` 与 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import read_phd_sequence

result = read_phd_sequence("read.phd", verbose=False)
print(result["file_name"], result["sequence"], result["qualities"], result["peak_locations"])
```

- **示例输出：**

```text
read1 ACGT [10, 20, 30, 40] [1, 4, 2, 3]
```

## 13) `dnakit.sequence.read_genepop_genotypes()`

- **作用：** 读取 GenePop 基因型文件，得到位点、群体和每个个体的两等位基因基因型，用于群体遗传统计。
- **处理方法：** 调用 `Bio.PopGen.GenePop.read`。缺失基因型 `0000` 返回为 `[null, null]`，其他基因型保持文件中的等位基因编码。
- **API：** `dnakit.sequence.read_genepop_genotypes(path[required], verbose[optional])`。
- **输入：** GenePop 文件路径（每个个体为两个等位基因）。
- **输出：** `loci`、`marker_length`、`population_count`、`individual_count`、`populations`（含个体名与基因型）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import read_genepop_genotypes

result = read_genepop_genotypes("pop.genepop", verbose=False)
print(result["loci"], result["population_count"], result["individual_count"])
print(result["populations"][0]["individuals"][1])
```

- **示例输出：**

```text
['Locus1', 'Locus2'] 2 4
{'name': 'ind2', 'genotypes': [[1, 2], [None, None]]}
```

## 14) `dnakit.sequence.read_gaf_annotations()`

- **作用：** 读取 GAF 1.0 或 2.x 基因注释文件，得到每条注释记录，并统计 GO 术语与 aspect。
- **处理方法：** 调用 `Bio.UniProt.GOA.gafiterator`。空的限定词、参考文献和同义词会被删除；记录数超过 `max_records` 时报错，而不是截断。
- **API：** `dnakit.sequence.read_gaf_annotations(path[required], max_records[optional], verbose[optional])`。
- **输入：** GAF 文件路径；`max_records` 默认 100,000。
- **输出：** `record_count`、`go_term_counts`、`aspect_counts`、`records` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import read_gaf_annotations

result = read_gaf_annotations("annotations.gaf", verbose=False)
print(result["record_count"], result["go_term_counts"], result["aspect_counts"])
```

- **示例输出：**

```text
2 {'GO:0003700': 1, 'GO:0006355': 1} {'F': 1, 'P': 1}
```

## 15) `dnakit.sequence.parse_kegg_gene_entries()`

- **作用：** 解析已保存的 KEGG GENE 平面文件，得到基因名、直系同源、物种、染色体位置和数据库交叉引用。
- **处理方法：** 调用 `Bio.KEGG.Gene.parse`。本函数只读取本地文件，不访问 KEGG REST 服务。
- **API：** `dnakit.sequence.parse_kegg_gene_entries(path[required], verbose[optional])`。
- **输入：** KEGG GENE 平面文件路径。
- **输出：** `entry_count`、`entries`（每条含 `entry`、`names`、`orthology`、`organism`、`position`、`dblinks`）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import parse_kegg_gene_entries

result = parse_kegg_gene_entries("gene.txt", verbose=False)
entry = result["entries"][0]
print(entry["entry"], entry["organism"]["code"], entry["orthology"])
```

- **示例输出：**

```text
hsa:7157 hsa [{'id': 'K04442', 'name': 'tumor protein p53'}]
```

## 16) `dnakit.sequence.load_substitution_matrix()`

- **作用：** 加载 Biopython 自带的替换矩阵，包括 NUC.4.4 等 DNA 比对打分矩阵。
- **处理方法：** 调用 `Bio.Align.substitution_matrices.load`。名称必须在 Biopython 的列表中；行与列按矩阵字母表的顺序输出。
- **API：** `dnakit.sequence.load_substitution_matrix(name[required], verbose[optional])`。
- **输入：** 矩阵名称，例如 `NUC.4.4` 或 `BLOSUM62`。
- **输出：** `name`、`alphabet`、`size`、`symmetric`、`matrix` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import load_substitution_matrix

result = load_substitution_matrix("NUC.4.4", verbose=False)
print(result["alphabet"], result["symmetric"])
print(result["matrix"][0][:4])
```

- **示例输出：**

```text
ATGCSWRYKMBVHDN True
[5.0, -4.0, -4.0, -4.0]
```

## 17) `dnakit.sequence.convert_quality_score()`

- **作用：** 在 Solexa 质量值与 PHRED 质量值之间换算单个分数，用于处理旧版 FASTQ 数据。
- **处理方法：** 调用 `Bio.SeqIO.QualityIO` 中的公式。Phred = 10·log10(10^(Q/10) + 1)；Solexa 换算使用相反的公式。PHRED 输入不能为负。
- **API：** `dnakit.sequence.convert_quality_score(value[required], direction[optional], verbose[optional])`。
- **输入：** 一个有限数值；`direction` 为 `solexa_to_phred`（默认）或 `phred_to_solexa`。
- **输出：** `direction`、`input`、`output` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import convert_quality_score

print(round(convert_quality_score(10.0, verbose=False)["output"], 4))
print(round(convert_quality_score(10.0, direction="phred_to_solexa", verbose=False)["output"], 4))
```

- **示例输出：**

```text
10.4139
9.5424
```

## 18) `dnakit.sequence.read_psl_alignments()`

- **作用：** 读取 PSL 格式的基因组两两比对（例如 BLAT 或 minimap2 的输出），得到目标与查询序列上的区块坐标。
- **处理方法：** 调用 `Bio.Align.parse(format='psl')`。`target_coordinates` 与 `query_coordinates` 为各序列上的区块边界；相邻坐标差相同的区块无空位。
- **API：** `dnakit.sequence.read_psl_alignments(path[required], max_alignments[optional], verbose[optional])`。
- **输入：** PSL 文件路径；`max_alignments` 默认 100,000。
- **输出：** `alignment_count`、`alignments`（每条含 `target_id`、`query_id`、起止坐标和区块坐标）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import read_psl_alignments

result = read_psl_alignments("alignments.psl", verbose=False)
print(result["alignment_count"])
print(result["alignments"][0])
```

- **示例输出：**

```text
1
{'target_id': 'chr1', 'query_id': 'read1', 'target_start': 100, 'target_end': 120, 'query_start': 0, 'query_end': 20, 'target_coordinates': [100, 120], 'query_coordinates': [0, 20]}
```

## 19) `dnakit.sequence.bootstrap_consensus_tree()`

- **作用：** 对一组比对好的 DNA 序列做自助法重采样，并返回多数规则一致树（Newick 格式）。
- **处理方法：** 调用 `Bio.Phylo.Consensus.bootstrap_consensus`，建树方法为 NJ 或 UPGMA，距离模型为 identity、blastn 或 trans。调用期间固定随机种子，结束后恢复原随机状态。
- **API：** `dnakit.sequence.bootstrap_consensus_tree(sequences[required], replicates[optional], method[optional], distance_model[optional], seed[optional], max_length[optional], verbose[optional])`。
- **输入：** ID 到等长 ACGT 比对序列的映射，至少三条，最多 200 条；`replicates` 为 1 到 1,000。
- **输出：** `ids`、`method`、`distance_model`、`replicates`、`seed`、`newick` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import bootstrap_consensus_tree

rows = {
    "s1": "ACGTACGTAAGGCCTTACGT",
    "s2": "ACGTACCTAAGGCATTACGA",
    "s3": "ACGTTCGTAAGGCCTAACGT",
    "s4": "ACGTTCGAAAGGCCTAACGA",
}
result = bootstrap_consensus_tree(rows, replicates=5, seed=1, verbose=False)
print(result["newick"])
```

- **示例输出：**

```text
(s3:0.025,s4:0.075,(s1:0.025,s2:0.125)100.00:0.08):0;
```

## 20) `dnakit.sequence.parse_primer3_output()`

- **作用：** 读取一份 Primer3 引物设计报告，得到引物对、产物长度和可选的内部探针。
- **处理方法：** 调用 `Bio.Emboss.Primer3.read`。只解析已保存的报告，不运行 Primer3；报告中只能包含一个目标序列。
- **API：** `dnakit.sequence.parse_primer3_output(path[required], verbose[optional])`。
- **输入：** Primer3 输出文件路径。
- **输出：** `pair_count`、`primer_pairs`（每对含产物长度、正向与反向引物、可选的 `internal`）、`comments` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import parse_primer3_output

result = parse_primer3_output("primer3.txt", verbose=False)
pair = result["primer_pairs"][0]
print(pair["product_size"], pair["forward"]["sequence"], pair["reverse"]["start"])
```

- **示例输出：**

```text
60 ACGTACGTACGTACGTACGT 60
```

## 21) `dnakit.sequence.nucleotide_distance()`

- **作用：** 计算两条等长、无空位的比对 DNA 序列之间的进化距离，可选 p 距离、JC69、K2P、F81、F84、TN93、LogDet 或 paralinear。
- **处理方法：** 调用 `skbio.sequence.distance` 中对应的函数。若 scikit-bio 返回 NaN（例如距离饱和），本函数返回 `distance: null` 且 `status: undefined`。
- **API：** `dnakit.sequence.nucleotide_distance(seq1[required], seq2[required], model[optional], verbose[optional])`。
- **输入：** 两条等长的 IUPAC DNA 字符串，不允许空位符号 `-`；`model` 为指定的模型名。
- **输出：** `model`、`sequence_length`、`distance`、`status` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import nucleotide_distance

left = "ACGTACGTAAGGCCTTACGT"
right = "ACGTACCTAAGGCATTACGA"
for model in ["pdist", "jc69", "tn93"]:
    result = nucleotide_distance(left, right, model=model, verbose=False)
    print(model, round(result["distance"], 6))
```

- **示例输出：**

```text
pdist 0.15
jc69 0.167358
tn93 0.170588
```

## 22) `dnakit.sequence.nucleotide_transition_matrix()`

- **作用：** 由一个进化距离得到 4 x 4 的核苷酸替换概率矩阵，用于模拟或似然计算。
- **处理方法：** 调用 `skbio.sequence.transition` 中的 JC69、K2P、F81、HKY85 或 TN93 函数，并使用 float64 精度。转换参数 `kappa` 必须在 (0, 1] 内，这是 scikit-bio 的约定。
- **API：** `dnakit.sequence.nucleotide_transition_matrix(distance[required], model[optional], kappa[optional], kappa_purine[optional], kappa_pyrimidine[optional], base_frequencies[optional], verbose[optional])`。
- **输入：** 非负距离；`base_frequencies` 为 A、C、G、T 的四个正数，和为 1（默认均为 0.25）。
- **输出：** `model`、`distance`、`alphabet`（`ACGT`）、`base_frequencies`、`matrix` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import nucleotide_transition_matrix

result = nucleotide_transition_matrix(0.2, verbose=False)
print(result["alphabet"], [round(value, 4) for value in result["matrix"][0]])
```

- **示例输出：**

```text
ACGT [0.8244, 0.0585, 0.0585, 0.0585]
```

## 23) `dnakit.sequence.alignment_distance_matrix()`

- **作用：** 计算 DNA 多序列比对各行两两之间的距离矩阵。
- **处理方法：** 调用 `skbio.alignment.align_dists`。输入的空位符号 `-` 会保留，并由 scikit-bio 的 `shared_by_all` 默认规则处理。
- **API：** `dnakit.sequence.alignment_distance_matrix(aligned_sequences[required], metric[optional], verbose[optional])`。
- **输入：** ID 到等长比对序列的映射（允许 `-`），至少两条，最多 500 条；`metric` 为 hamming、pdist、jc69、k2p、f81、f84、tn93、logdet 或 paralin。
- **输出：** `metric`、`ids`、`matrix`（对称方阵）和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import alignment_distance_matrix

rows = {"x": "ACGTACGTAA", "y": "ACGTTCGTAC", "z": "ACGAACGTAA"}
result = alignment_distance_matrix(rows, metric="hamming", verbose=False)
print(result["ids"])
print(result["matrix"])
```

- **示例输出：**

```text
['x', 'y', 'z']
[[0.0, 0.2, 0.1], [0.2, 0.0, 0.3], [0.1, 0.3, 0.0]]
```

## 24) `dnakit.sequence.multiple_nucleotide_alignment()`

- **作用：** 对数条 DNA 序列做渐进式多序列比对，不需要 MAFFT、MUSCLE 或 Clustal 等外部程序。
- **处理方法：** 调用 `skbio.alignment.multi_align_nucl`，匹配、错配与空位罚分可以设置。输出按输入 ID 的字典序排列；去掉空位后应与输入序列一致。
- **API：** `dnakit.sequence.multiple_nucleotide_alignment(sequences[required], match_score[optional], mismatch_score[optional], gap_cost[optional], free_ends[optional], verbose[optional])`。
- **输入：** ID 到未比对 IUPAC DNA 序列的映射，2 到 100 条，每条长度不超过 100,000。
- **输出：** `ids`、`alignment`（ID 到带空位序列）、`alignment_length`、罚分参数和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import multiple_nucleotide_alignment

rows = {"s1": "ACGTACGTAAGGCC", "s2": "ACGTACCTAAGGCA", "s3": "ACGTTCGTAAGGCCTA"}
result = multiple_nucleotide_alignment(rows, verbose=False)
print(result["alignment_length"])
for identifier, row in result["alignment"].items():
    print(identifier, row)
```

- **示例输出：**

```text
16
s1 ACGTACGTAAGGCC--
s2 ACGTACCTAAGGCA--
s3 ACGTTCGTAAGGCCTA
```

## 25) `dnakit.sequence.alignment_score()`

- **作用：** 按匹配、错配和空位罚分计算一组 DNA 比对的得分，用于比较不同比对的质量。
- **处理方法：** 调用 `skbio.alignment.align_score`，输入为 ID 到等长比对序列的映射。`free_ends` 控制两端空位是否计罚。
- **API：** `dnakit.sequence.alignment_score(aligned_sequences[required], match_score[optional], mismatch_score[optional], gap_cost[optional], free_ends[optional], verbose[optional])`。
- **输入：** ID 到等长比对序列的映射（允许 `-`），至少两条。
- **输出：** `ids`、`alignment_length`、`score`、罚分参数和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import alignment_score

print(alignment_score({"x": "ACGTA", "y": "ACGTA"}, verbose=False)["score"])
print(alignment_score({"x": "AC-TA", "y": "ACGTA"}, verbose=False)["score"])
```

- **示例输出：**

```text
5.0
2.0
```

## 26) `dnakit.sequence.sniff_sequence_format()`

- **作用：** 检测本地序列文件的格式，并给出 scikit-bio 建议的读取选项，用于在读取前确认文件类型。
- **处理方法：** 调用 `skbio.io.sniff`，只读取文件开头部分，不转换文件。无法识别时返回配置错误。
- **API：** `dnakit.sequence.sniff_sequence_format(path[required], verbose[optional])`。
- **输入：** 本地序列文件路径。
- **输出：** `format`、`reader_options`、`source` 和 `provenance`。
- **示例代码：**

```python
from dnakit.sequence import sniff_sequence_format

result = sniff_sequence_format("reads.fq", verbose=False)
print(result["format"], result["reader_options"])
```

- **示例输出：**

```text
fastq {}
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-formats]"   # Biopython
python -m pip install "dnakit[bio-skibio]"    # scikit-bio
```

## 验证方式

单元测试位于 `others/tests/unit/test_full_G1a.py`，每个函数至少一个测试，使用本地临时文件并调用真实的上游函数。
