# Biopython 系统发育、NEXUS、树查询、dN/dS 与 PAML

这一组 `dnakit.phylogeny` 函数封装了 [Biopython](https://biopython.org)（`Bio.Nexus`、`Bio.Phylo`、`Bio.Align.analysis`、`Bio.codonalign`、`Bio.Phylo.PAML`；测试版本 Biopython 1.88）和 PAML 程序。比对是“分类单元名到比对序列”的映射；树是 Newick 字符串。Biopython 只在函数运行时导入，可用 `python -m pip install "dnakit[bio-formats]"` 安装。含随机性的函数接受 `seed`，调用结束后恢复全局随机状态。`run_paml` 需要 PAML（`codeml`、`baseml`、`yn00`）在 `PATH` 中，例如 `conda install -c bioconda paml`。示例从模块 `dnakit.phylogeny.biopython_p19_phylo_functions` 导入。

## 示例数据

`codons.phy` 是 PHYLIP 格式、4 条序列、600 个碱基、符合读码框的编码比对（模拟数据：一条随机序列和三条突变拷贝）。`pair_yn00.out` 是 PAML 4.10.10 对该比对中序列对 `a`、`d` 运行 `yn00` 得到的真实结果文件。

**示例文件：** `codons.phy`

```text
 4 600
a  GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGGCTGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCGCTTCGACTATGCTACTGCGGTTAACGGGGATGGCAAGTACATTTTTTCGTAGATGTGCCTTGCTAACGAAAGTATTAAACACGTCCCTCACAATAGAATCAGCTTTGGACGCGCGACGGCCGTTCCAGAAAATCTTTGAATACTCAATCCTGCGGGTTCGGGCTCCTAAAACCCATGCTTTGTGTTACCCAGTTCGAGCGCATAGGGAATTCAGGTCCACACATGGCTGGATCCCCATGATATTCAAGAACTATACATTAAGTGCTACCTCCAGAACACATGTTTCAGTCACGGCTTGCCATCATCGATCACGGAATGTAGCATCA
b  GGATCACAGTCTACACTGCTCACTCCAAACCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGAATGTATACCACTGGGGCTGATACGGCGGAGGCCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTCCATCTGCATGGAGAGGGTGGGCATGGATGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGTTTATCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCGCTTCGACTATGCTACTGCGGTTAACGGGGATGGCAAGTAAATTTTTTCGTGGATCTGCCTTGCTAACGAAAGTATTAAACACGTCCCTCACCGTAGAATCAGCTTTGGACGCGCGACGGCCGTTCCAGAAAATCTTTGAATACTCAATCCTGCGGGTTCGGGCACCTAAATCCCATGCTTTGTGTTACCCAGTTCGAGCGCATAGGGAATTCAGGTCCACACATGGCTGGATCCCGATGATATTCAAGAACTATACATTAAGTGCTACCGCCAGAACACATGTTTCAGTCACGGCTTGCCATCATCGATCACGGAATGTAGCATCA
c  GGATCACAGTTTACACTGCTGACTCCAACCCCGACCCCTGAGTCCGAGGAGAGGGTGCTTCAGCGTATGTATACCACTGGGGCTGATACAACGGAGAGCACGTCAAAACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCAAGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATAGTACCGATTGTAGAGAGGGGCTTGTCCTTCCCGAGAGCGTTTCTGTTTCGGTGTAGGTGCGCTTCGACAACGCTAATCCGGTTAACGGGGATGGCAACTACTTTTTATCGTAGATGTGCCTTGCTAACGAAAGTATTAAACACGTCCGTCACAATAGAATCAGCTTTGGACGCGCGATGGCCGTTCCAGGAAACCTTTGAATACTCAATCCTGCGGCTTCGGGCTCCTAAATCCCATGCCTTGTGTTTCCCAGTTCGAGCGCATAGGGTATTCAGGTCCACACATGGCTGTATCCCCATGTTATTCAAGAAATATACATTAAGTGCTACCTCCAGAACACATGTTTCAGTCACGGCTTGCCATCATCGATCACGGAATGTAGCATCA
d  GGATCACAGTCTACACGGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGTGGGTGCTTCAAACTATGTATACCACTCGGGCTGATACGTCGGAGGGGACGTTAATACGGTTCCATGCCATAATGCATGCTCTTGTGGTTTATCTGCAAGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCTGTGATCTGGACCTCCAAGCCACAGCTCACTGTACAGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCATTGCTGTTTCGGTGTAGGTGCGCTTCGACTATGCTTCTGCGGTTACCGGGAATGGCACGTGCATTTTTTCGTAGATGCTTCCTGCTAACGAAATTATTTAACACGTCGCTCACAATAGAATCAGCTTTGGAGTCGCGACGGCCGTTCAAGGCTATCTTTGAATACTCAATCCTGTGGGCCCGGCCTCCTATACCCCATGCTTTGTGTTACCCAGTTCCAGCGCATAGCGAATTCAGGTACACCCATGGCGCTATCACCATGATATTCAAGAACTATACAATAAGTGCTATCTCCAGAACACATGTTTCAGCCACAGCTTGCCATCATCGATCACGGAATGTAGCATCA
```

**示例文件：** `pair_yn00.out`

```text
YN00         seq.phy

ns =   2	ls = 200

Codon position x base (3x4) table for each sequence.

a                   
position  1:    T:0.23000    C:0.25500    A:0.26000    G:0.25500
position  2:    T:0.28000    C:0.31000    A:0.19500    G:0.21500
position  3:    T:0.25500    C:0.20000    A:0.22500    G:0.32000

d                   
position  1:    T:0.25000    C:0.25500    A:0.25500    G:0.24000
position  2:    T:0.28000    C:0.32000    A:0.19500    G:0.20500
position  3:    T:0.25500    C:0.20000    A:0.23000    G:0.31500

Average
position  1:    T:0.24000    C:0.25500    A:0.25750    G:0.24750
position  2:    T:0.28000    C:0.31500    A:0.19500    G:0.21000
position  3:    T:0.25500    C:0.20000    A:0.22750    G:0.31750

Codon usage for each species
----------------------------------------------------------------------
Phe TTT   5   5 | Ser TCT   1   1 | Tyr TAT   2   3 | Cys TGT   3   2
    TTC   5   6 |     TCC   6   4 |     TAC   2   3 |     TGC   3   4
Leu TTA   3   3 |     TCA   7   6 | *** TAA   0   0 | *** TGA   0   0
    TTG   4   4 |     TCG   1   4 |     TAG   0   0 | Trp TGG   4   5
----------------------------------------------------------------------
Leu CTT   2   3 | Pro CCT   2   4 | His CAT  11   9 | Arg CGT   1   2
    CTC   3   3 |     CCC   2   1 |     CAC   0   0 |     CGC   0   0
    CTA   3   1 |     CCA   3   4 | Gln CAA   0   2 |     CGA   3   2
    CTG   6   5 |     CCG   3   3 |     CAG   5   4 |     CGG   7   8
----------------------------------------------------------------------
Ile ATT   1   0 | Thr ACT   3   5 | Asn AAT   2   1 | Ser AGT   4   2
    ATC   4   5 |     ACC   5   5 |     AAC   2   2 |     AGC   0   1
    ATA   4   7 |     ACA   6   5 | Lys AAA   3   1 | Arg AGA   2   2
Met ATG   4   5 |     ACG   6   4 |     AAG   1   3 |     AGG   5   3
----------------------------------------------------------------------
Val GTT   4   3 | Ala GCT   8   9 | Asp GAT   1   1 | Gly GGT   1   1
    GTC   1   0 |     GCC   3   4 |     GAC   1   0 |     GGC   3   2
    GTA   4   3 |     GCA   2   4 | Glu GAA   3   3 |     GGA   2   3
    GTG   3   3 |     GCG   4   1 |     GAG   6   7 |     GGG   5   4
----------------------------------------------------------------------


Sums
------------------------------------------------------
Phe TTT  10 | Ser TCT   2 | Tyr TAT   5 | Cys TGT   5
    TTC  11 |     TCC  10 |     TAC   5 |     TGC   7
Leu TTA   6 |     TCA  13 | *** TAA   0 | *** TGA   0
    TTG   8 |     TCG   5 |     TAG   0 | Trp TGG   9
------------------------------------------------------
Leu CTT   5 | Pro CCT   6 | His CAT  20 | Arg CGT   3
    CTC   6 |     CCC   3 |     CAC   0 |     CGC   0
    CTA   4 |     CCA   7 | Gln CAA   2 |     CGA   5
    CTG  11 |     CCG   6 |     CAG   9 |     CGG  15
------------------------------------------------------
Ile ATT   1 | Thr ACT   8 | Asn AAT   3 | Ser AGT   6
    ATC   9 |     ACC  10 |     AAC   4 |     AGC   1
    ATA  11 |     ACA  11 | Lys AAA   4 | Arg AGA   4
Met ATG   9 |     ACG  10 |     AAG   4 |     AGG   8
------------------------------------------------------
Val GTT   7 | Ala GCT  17 | Asp GAT   2 | Gly GGT   2
    GTC   1 |     GCC   7 |     GAC   1 |     GGC   5
    GTA   7 |     GCA   6 | Glu GAA   6 |     GGA   5
    GTG   6 |     GCG   5 |     GAG  13 |     GGG   9
------------------------------------------------------



(A) Nei-Gojobori (1986) method



Nei & Gojobori 1986. dN/dS (dN, dS)
(Note: This matrix is not used in later ML. analysis.
Use runmode = -2 for ML pairwise comparison.)

a                   
d                    1.1566 (0.1033 0.0894)


(B) Yang & Nielsen (2000) method

Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43

(equal weighting of pathways)

seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE

   2    1   163.2   436.8   0.2998  1.1043  1.2262 0.1052 +- 0.0164  0.0858 +- 0.0242


(C) LWL85, LPB93 & LWLm methods

Li W.-H., C.-I. Wu, Luo (1985) A new method for estimating synonymous and nonsynonymous rates of nucleotide substitutions considering the relative likelihood of nucleotide and codon changes. Mol. Biol. Evol. 2: 150-174.
Li W-H (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol. 36:96-99
Pamilo P, Bianchi NO (1993) Evolution of the Zfx and Zfy genes - rates and interdependence between the genes. Mol. Biol. Evol. 10:271-281
Yang Z (2006) Computational Molecular Evolution. Oxford University Press, Oxford. Eqs. 2.12 & 2.13

2 (d) vs. 1 (a)

L(i):      378.5     113.0     108.5  sum=    600.0
Ns(i):    9.5000    3.5000    5.0000  sum=  18.0000
Nv(i):   23.5000    9.5000    5.0000  sum=  38.0000
A(i):     0.0264    0.0329    0.0502
B(i):     0.0663    0.0920    0.0483
LWL85:  dS =  0.0986 dN =  0.1002 w = 1.0165 S =  146.2 N =  453.8
LWL85m: dS =  0.0868 dN =  0.1048 w = 1.2080 S =  166.1 N =  433.9 (rho = 0.510)
LPB93:  dS =  0.0897 dN =  0.0986 w = 1.0991
```

## 1) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_terminal_gaps_to_missing()`

- **作用：** 把每条比对序列首尾的缺口替换为缺失数据符号。
- **处理方法：** 把比对写成 NEXUS 数据块，用 `Bio.Nexus.Nexus.Nexus` 读入并调用 `terminal_gap_to_missing`。内部缺口保持不变；`skip_n=False` 时末端的 N 也按缺口处理。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_terminal_gaps_to_missing(alignment, *, missing='?', skip_n=True, datatype='dna')`。
- **输入：** 分类单元名到比对序列的映射（缺口 `-`、缺失 `?`，长度相等）、缺失符号和 `datatype`（`dna`、`rna` 或 `nucleotide`）。
- **输出：** `taxa` 与转换后的 `alignment`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_terminal_gaps_to_missing

result = nexus_terminal_gaps_to_missing({"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"})
for name, row in result["alignment"].items():
    print(name, row)
```

- **示例输出：**

```text
A ??ACGT--A
B ?ACGTTTAA
C ?ACG--TTA
```

## 2) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_combine_matrices()`

- **作用：** 把多个字符矩阵串联成一个超矩阵。
- **处理方法：** 把每个比对读成 NEXUS 矩阵并调用 `Bio.Nexus.Nexus.combine`。某数据集中缺少的分类单元在其列上填 `?`；每个数据集生成一个记录其列的字符集。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_combine_matrices(matrices, *, datatype='dna')`。
- **输入：** 数据集名到比对（分类单元名到比对序列）的映射，以及 `datatype`。
- **输出：** `taxa`、`character_count`、合并后的 `alignment`、`charsets` 和 `partition`（0 起始的列号）。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_combine_matrices

result = nexus_combine_matrices({
    "g1": {"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"},
    "g2": {"A": "ACGTAA", "D": "TTTTGG"},
})
print(result["taxa"], result["character_count"])
for name, row in result["alignment"].items():
    print(name, row)
print(result["charsets"])
```

- **示例输出：**

```text
['A', 'B', 'C', 'D'] 15
A --ACGT--AACGTAA
B ?ACGTTTAA??????
C -ACG--TTA??????
D ?????????TTTTGG
{'g1': [0, 1, 2, 3, 4, 5, 6, 7, 8], 'g2': [9, 10, 11, 12, 13, 14]}
```

## 3) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_write_partition_files()`

- **作用：** 为比对的每个字符分区写出一个 NEXUS 数据文件。
- **处理方法：** 用分区（例如密码子位置）调用 `Nexus.write_nexus_data_partitions`。每个文件 `<file_stem>_<partition>.data` 只含该分区的列；`exclude` 去掉列，`delete` 去掉分类单元。输出目录必须尚不存在。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_write_partition_files(alignment, charpartition, output_dir, *, file_stem='alignment', exclude=None, delete=None, interleave=False, datatype='dna')`。
- **输入：** 比对映射、`charpartition`（分区名到 0 起始列号）、新的 `output_dir`、`file_stem`、`exclude`、`delete`、`interleave` 和 `datatype`。
- **输出：** `files`（分区名到路径）与 `partition_sizes`。
- **示例代码：**

```python
import os, tempfile
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_write_partition_files

target = os.path.join(tempfile.mkdtemp(), "parts")
result = nexus_write_partition_files(
    {"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"},
    {"first": [0, 1, 2], "rest": [3, 4, 5, 6, 7, 8]},
    target,
)
print(result["partition_sizes"])
print(sorted(os.path.basename(path) for path in result["files"].values()))
print(open(result["files"]["first"]).read())
```

- **示例输出：**

```text
{'first': 3, 'rest': 6}
['alignment_first.data', 'alignment_rest.data']
#NEXUS
[
Partition: first
]
begin data;
dimensions ntax=3 nchar=3;
format datatype=dna missing=? gap=-;
matrix
A --A
B ?AC
C -AC
;
end;
```

## 4) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_weighted_stepmatrix()`

- **作用：** 根据比对中观察到的字符状态计算加权简约步长矩阵。
- **处理方法：** 遵循 Wheeler (1990)：统计两个状态在同一位点共现的次数，归一化后取负对数（`Bio.Nexus.Nexus.Nexus.weighted_stepmatrix`）。从未共现的状态对记为 `inf`。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_weighted_stepmatrix(alignment, *, name='your_name_here', exclude=None, delete=None, datatype='dna')`。
- **输入：** 至少含两个分类单元的比对映射、矩阵名称 `name`、`exclude`（0 起始列号）、`delete`（分类单元）和 `datatype`。
- **输出：** `symbols`、NEXUS 格式的 `stepmatrix` 文本和 `weights`（`null` 表示 `inf`）。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_weighted_stepmatrix

result = nexus_weighted_stepmatrix({"A": "--ACGT--A", "B": "?ACGTTTAA", "C": "-ACG--TTA"}, name="steps")
print(result["symbols"])
print("\n".join(line.rstrip() for line in result["stepmatrix"].splitlines()))
```

- **示例输出：**

```text
['A', 'C', 'G', 'T', '-']
usertype steps stepmatrix=5
        A        C        G        T        -
[A]     .       1.39     inf.     1.39     inf.
[C]    1.39      .       1.39     inf.     inf.
[G]    inf.     1.39      .       1.39     inf.
[T]    1.39     inf.     1.39      .       inf.
[-]    inf.     inf.     inf.     inf.      .
;
```

## 5) `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_bootstrap_alignment()`

- **作用：** 对比对的列做有放回重抽样。
- **处理方法：** 去掉排除的列和删除的分类单元后调用 `Bio.Nexus.Nexus.Nexus.bootstrap`。给定 `seed` 时结果可复现；调用结束后恢复全局随机状态。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.nexus_bootstrap_alignment(alignment, *, seed=None, exclude=None, delete=None, datatype='dna')`。
- **输入：** 比对映射、`seed`、`exclude`、`delete` 和 `datatype`。
- **输出：** 自举后的 `alignment`（长度与去除排除列后的输入相同）。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import nexus_bootstrap_alignment

result = nexus_bootstrap_alignment({"A": "ACGTACGTAC", "B": "ACGTACGTTC", "C": "ACGAACGTTC"}, seed=7)
for name, row in result["alignment"].items():
    print(name, row)
```

- **示例输出：**

```text
A CGGACACCCA
B CGGACTCCCA
C CGGACTCCCA
```

## 6) `dnakit.phylogeny.biopython_p19_phylo_functions.biopython_alignment_distance_matrix()`

- **作用：** 用 Biopython 的 `DistanceCalculator` 计算 DNA 比对的两两距离矩阵。
- **处理方法：** `identity` 为比对位点中不同位点的比例；`blastn`、`megablast`、`trans` 和 `hoxd70` 由核苷酸替换矩阵导出距离（`Bio.Phylo.TreeConstruction.DistanceCalculator.get_distance`）。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.biopython_alignment_distance_matrix(alignment, *, model='identity')`。
- **输入：** 至少含两条序列的比对映射和 `model`。
- **输出：** `ids` 和完整的对称 `matrix`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import biopython_alignment_distance_matrix

aln = {"a": "ACGTACGTAC", "b": "ACGTACGTTC", "c": "ACGAACGTTC", "d": "TCGAACGATC"}
for model in ("identity", "blastn"):
    result = biopython_alignment_distance_matrix(aln, model=model)
    print(model, result["ids"], [round(value, 3) for value in result["matrix"][3]])
```

- **示例输出：**

```text
identity ['a', 'b', 'c', 'd'] [0.4, 0.3, 0.2, 0.0]
blastn ['a', 'b', 'c', 'd'] [0.72, 0.54, 0.36, 0.0]
```

## 7) `dnakit.phylogeny.biopython_p19_phylo_functions.parsimony_tree_search()`

- **作用：** 用最近邻交换搜索最大简约树。
- **处理方法：** 运行 Biopython 的 `ParsimonyTreeConstructor`，搜索器为 `NNITreeSearcher`，评分器为 `ParsimonyScorer`（Fitch，单位代价）。没有起始树时从基于一致性距离的邻接树开始。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.parsimony_tree_search(alignment, *, starting_tree=None)`。
- **输入：** 至少含三条序列的比对映射，以及可选的、带分支长度且分类单元相同的 Newick `starting_tree`。
- **输出：** `newick` 树及其 `parsimony_score`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import parsimony_tree_search

result = parsimony_tree_search({"a": "ACGTACGTAC", "b": "ACGTACGTTC", "c": "ACGAACGTTC", "d": "TCGAACGATC"})
print(result["newick"])
print(result["parsimony_score"], result["starting_tree"])
```

- **示例输出：**

```text
(d:0.1625,((a:0.075,b:0.05)Inner1:0.025,c:0.05)Inner2:0.0875)Inner3:0;
4.0 neighbour joining (identity distance)
```

## 8) `dnakit.phylogeny.biopython_p19_phylo_functions.strict_consensus_tree()`

- **作用：** 构建多棵树的严格一致树：只保留所有树中都出现的分支。
- **处理方法：** 按末端名称匹配分支。Biopython 的 `strict_consensus` 比较末端位置，只有各树末端顺序一致时才正确；本函数在这种情况下得到相同的树，其他情况下也正确。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.strict_consensus_tree(newick_trees)`。
- **输入：** 至少两个分类单元相同的 Newick 字符串列表。
- **输出：** `newick`、`tree_count` 和 `shared_clade_count`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import strict_consensus_tree

result = strict_consensus_tree(["((A,B),(C,D),E);", "((A,B),(C,E),D);"])
print(result["newick"], result["shared_clade_count"])
```

- **示例输出：**

```text
((A:0,B:0):0,C:0,D:0,E:0):0; 1
```

## 9) `dnakit.phylogeny.biopython_p19_phylo_functions.collapse_tree_genera()`

- **作用：** 把末端同属一个属的子树折叠成一个末端。
- **处理方法：** 属名取分类单元名的第一个词（`Genus_species`）。使用 `Bio.Nexus.Trees.Tree.collapse_genera`；折叠后的子树成为名为 `<Genus> <collapsed>` 的单个末端。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.collapse_tree_genera(newick, *, space_equals_underscore=True)`。
- **输入：** Newick 字符串和 `space_equals_underscore`。
- **输出：** 折叠后的 `newick` 以及折叠前后的末端数。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import collapse_tree_genera

result = collapse_tree_genera("((Homo_sapiens:1,Homo_erectus:1):1,(Pan_troglodytes:1,Gorilla_gorilla:1):1);")
print(result["newick"])
print(result["terminal_count_before"], result["terminal_count_after"])
```

- **示例输出：**

```text
(Homo <collapsed>:1.00000,(Pan_troglodytes:1.00000,Gorilla_gorilla:1.00000)0.00000:1.00000);
4 3
```

## 10) `dnakit.phylogeny.biopython_p19_phylo_functions.split_tree_clade()`

- **作用：** 把指定名称的分支拆成新的子分支。
- **处理方法：** 调用 `Bio.Phylo.BaseTree.TreeMixin.split`。新子分支命名为 `<target>0`、`<target>1` ……并获得 `branch_length`；末端会变成内部节点。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.split_tree_clade(newick, target, *, count=2, branch_length=1.0)`。
- **输入：** Newick 字符串、唯一的 `target` 分支名、`count` 和 `branch_length`。
- **输出：** `new_clades` 和编辑后的 `newick`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import split_tree_clade

result = split_tree_clade("((A:1,B:1):1,C:2);", "C", count=3, branch_length=0.5)
print(result["new_clades"])
print(result["newick"])
```

- **示例输出：**

```text
['C0', 'C1', 'C2']
((A:1,B:1):1,(C0:0.5,C1:0.5,C2:0.5)C:2):0;
```

## 11) `dnakit.phylogeny.biopython_p19_phylo_functions.random_tree()`

- **作用：** 在给定分类单元上生成随机二叉树。
- **处理方法：** 调用 `Bio.Phylo.BaseTree.Tree.randomized`：反复拆分随机末端并打乱名称。它不是溯祖模型或生灭模型。给定 `seed` 时树可复现。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.random_tree(taxa, *, branch_length=1.0, branch_stdev=None, seed=None)`。
- **输入：** 2 到 10000 个互不相同的分类单元名、`branch_length`、可选的 `branch_stdev` 和 `seed`。
- **输出：** `newick` 树。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import random_tree

print(random_tree(["A", "B", "C", "D", "E"], seed=3)["newick"])
print(random_tree(["A", "B", "C", "D", "E"], seed=3)["newick"])
```

- **示例输出：**

```text
(((C:1,(A:1,E:1)n001:1)n00:1,D:1)n0:1,B:1):0;
(((C:1,(A:1,E:1)n001:1)n00:1,D:1)n0:1,B:1):0;
```

## 12) `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_bifurcating()`

- **作用：** 检验一棵树是否严格二叉。
- **处理方法：** 调用 `Bio.Phylo.BaseTree.TreeMixin.is_bifurcating`：每个内部节点有两个子节点；根有三个子节点（无根形式）也被接受。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_bifurcating(newick)`。
- **输入：** Newick 字符串。
- **输出：** `is_bifurcating`、`terminal_count` 以及内部节点子节点数的集合。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import tree_is_bifurcating

for text in ("((A,B),(C,D));", "(A,B,C,D);", "((A,B,C),D);"):
    result = tree_is_bifurcating(text)
    print(text, result["is_bifurcating"], result["internal_child_counts"])
```

- **示例输出：**

```text
((A,B),(C,D)); True [2]
(A,B,C,D); False [4]
((A,B,C),D); False [2, 3]
```

## 13) `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_parent_of()`

- **作用：** 检验一个命名分支是否包含另一个命名分支。
- **处理方法：** 对分支 `ancestor` 以分支 `descendant` 调用 `Bio.Phylo.BaseTree.TreeMixin.is_parent_of`。与 Biopython 一致，分支被视为自己的后代。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.tree_is_parent_of(newick, ancestor, descendant)`。
- **输入：** Newick 字符串和两个唯一的分支名。
- **输出：** `is_parent_of`（真或假）。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import tree_is_parent_of

tree = "((A,B)AB,(C,D)CD)root;"
print(tree_is_parent_of(tree, "AB", "A")["is_parent_of"])
print(tree_is_parent_of(tree, "AB", "C")["is_parent_of"])
```

- **示例输出：**

```text
True
False
```

## 14) `dnakit.phylogeny.biopython_p19_phylo_functions.find_tree_clades()`

- **作用：** 按名称、是否为末端或分支长度查找分支。
- **处理方法：** 使用 `Bio.Phylo` 的 `find_clades`（全部匹配）或 `find_any`（`first_only=True`）。`name` 是必须匹配整个名称的正则表达式；`terminal` 选择末端或内部节点；`branch_length` 按相等比较（容差 1e-9）。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.find_tree_clades(newick, *, name=None, terminal=None, branch_length=None, order='preorder', first_only=False, max_results=10000)`。
- **输入：** Newick 字符串、`name`、`terminal`、`branch_length`、`order`（`preorder`、`postorder` 或 `level`）、`first_only` 和 `max_results`。
- **输出：** `match_count` 和 `clades` 列表（名称、分支长度、置信度、是否末端、末端数）。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import find_tree_clades

tree = "((A_x:1,A_y:2)0.9:1,(C:1,D:1)0.8:2,E:3);"
print([c["name"] for c in find_tree_clades(tree, name="A.*")["clades"]])
print(find_tree_clades(tree, terminal=False)["match_count"])
print(find_tree_clades(tree, branch_length=3)["clades"])
```

- **示例输出：**

```text
['A_x', 'A_y']
3
[{'name': 'E', 'branch_length': 3.0, 'confidence': None, 'terminal': True, 'terminal_count': 1}]
```

## 15) `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_identical()`

- **作用：** 检验两棵树的拓扑和分类单元是否相同。
- **处理方法：** 调用 `Bio.Nexus.Trees.Tree.is_identical`，比较子节点的嵌套集合；忽略分支长度和子节点顺序。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_identical(newick_a, newick_b)`。
- **输入：** 两个 Newick 字符串。
- **输出：** `identical`（真或假）。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import trees_are_identical

print(trees_are_identical("((A,B),(C,D));", "((D,C),(B,A));")["identical"])
print(trees_are_identical("((A,B),(C,D));", "((A,C),(B,D));")["identical"])
```

- **示例输出：**

```text
True
False
```

## 16) `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_compatible()`

- **作用：** 检查两棵树中是否有相互冲突的受支持分支。
- **处理方法：** 调用 `Bio.Nexus.Trees.Tree.is_compatible`：比较支持度不低于 `threshold` 的内部分支，相互重叠又不嵌套的分支对即为冲突。支持度须写成 `)支持度:长度`，例如 `(A,B)90:1`。`strict=True` 时分类单元必须相同。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.trees_are_compatible(newick_a, newick_b, *, threshold, strict=True)`。
- **输入：** 两个 Newick 字符串、`threshold` 和 `strict`。
- **输出：** `compatible`、`conflict_count` 和 `conflicts`（分支、支持度与共有分类单元）。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import trees_are_compatible

same = trees_are_compatible("((A,B)90:1,(C,D)90:1,E);", "((A,B)95:1,(C,D)80:1,E);", threshold=50)
print(same["compatible"], same["conflict_count"])
clash = trees_are_compatible("((A,B)90:1,(C,D)90:1,E);", "((A,C)90:1,(B,D)90:1,E);", threshold=50)
print(clash["compatible"], clash["conflict_count"])
print(clash["conflicts"][0])
```

- **示例输出：**

```text
True 0
False 4
{'clade_a': ['A', 'B'], 'support_a': 90.0, 'clade_b': ['A', 'C'], 'support_b': 90.0, 'shared': ['A']}
```

## 17) `dnakit.phylogeny.biopython_p19_phylo_functions.convert_tree_support()`

- **作用：** 转换分支支持度。
- **处理方法：** `absolute_to_relative` 把分支计数（如 PHYLIP `consense` 的输出）除以重复次数（`convert_absolute_support`）。`branchlength_to_support` 把被读成分支长度的数值（例如 PAUP 自举树）移入支持度字段（`branchlength2support`）。两者都是 `Bio.Nexus.Trees.Tree` 的方法。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.convert_tree_support(newick, *, mode='absolute_to_relative', replicates=None)`。
- **输入：** Newick 字符串、`mode` 和 `replicates`（`absolute_to_relative` 必需）。
- **输出：** 带 `支持度:长度` 标签的 `newick`，以及含分类单元、支持度和分支长度的 `nodes`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import convert_tree_support

result = convert_tree_support("((A:1,B:2)90:1,(C:1,D:1)80:2,E:3);", replicates=100)
print(result["newick"])
moved = convert_tree_support("((A:0.9,B:2)0.9:1,(C:1,D:1)0.8:2,E:3);", mode="branchlength_to_support")
print(moved["newick"])
```

- **示例输出：**

```text
((A:1.00000,B:2.00000)0.90:1.00000,(C:1.00000,D:1.00000)0.80:2.00000,E:3.00000);
((A:0.00000,B:0.00000)1.00:0.00000,(C:0.00000,D:0.00000)2.00:0.00000,E:0.00000);
```

## 18) `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_estimates()`

- **作用：** 估计密码子比对的两两 dN 和 dS。
- **处理方法：** 调用 `Bio.Align.analysis.calculate_dn_ds_matrix`。`NG86` 为 Nei-Gojobori，`LWL85` 为 Li-Wu-Luo，`YN00` 为 Yang-Nielsen，`ML` 为 Goldman 与 Yang 的最大似然法。饱和的序列对（Biopython 返回 -1）报告为 `null`。默认的 `YN00` 与 `pairwise_dn_ds_matrix` 默认的 `NG86` 不同。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_estimates(codon_alignment, *, method='YN00', genetic_code=1)`。
- **输入：** 按读码框对齐、无缺口（仅 A、C、G、T）的编码序列比对映射、`method` 和 `genetic_code`（NCBI 密码子表编号）。
- **输出：** `ids` 以及完整对称的 `dn` 和 `ds` 矩阵。
- **示例代码：**

```python
rows = dict(line.split() for line in open("codons.phy").read().splitlines()[1:])
from dnakit.phylogeny.biopython_p19_phylo_functions import codon_dn_ds_estimates

for method in ("NG86", "LWL85", "YN00"):
    result = codon_dn_ds_estimates(rows, method=method)
    print(method, "dN(a,d) =", round(result["dn"][0][3], 4), "dS(a,d) =", round(result["ds"][0][3], 4))
```

- **示例输出：**

```text
NG86 dN(a,d) = 0.1017 dS(a,d) = 0.0936
LWL85 dN(a,d) = 0.0981 dS(a,d) = 0.1017
YN00 dN(a,d) = 0.0992 dS(a,d) = 0.1019
```

## 19) `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_trees()`

- **作用：** 由密码子比对构建 dN 树和 dS 树。
- **处理方法：** 调用 `Bio.codonalign.codonalignment.CodonAlignment.get_dn_ds_tree`：把两两 dN、dS 矩阵（`NG86`、`LWL85`、`YN00` 或 `ML`）用 UPGMA 或邻接法建成树。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.codon_dn_ds_trees(codon_alignment, *, dn_ds_method='NG86', tree_method='UPGMA', genetic_code=1)`。
- **输入：** 如上的密码子比对映射、`dn_ds_method`、`tree_method`（`UPGMA` 或 `NJ`）和 `genetic_code`。
- **输出：** `dn_tree_newick` 和 `ds_tree_newick`。
- **示例代码：**

```python
rows = dict(line.split() for line in open("codons.phy").read().splitlines()[1:])
from dnakit.phylogeny.biopython_p19_phylo_functions import codon_dn_ds_trees

result = codon_dn_ds_trees(rows, dn_ds_method="NG86", tree_method="UPGMA")
print(result["dn_tree_newick"])
print(result["ds_tree_newick"])
```

- **示例输出：**

```text
(d:0.066962759,(c:0.033571855,(b:0.010170115,a:0.010170115)Inner1:0.02340174)Inner2:0.033390904)Inner3:0;
(d:0.066604305,(b:0.03646134,(c:0.02201535,a:0.02201535)Inner1:0.01444599)Inner2:0.030142964)Inner3:0;
```

## 20) `dnakit.phylogeny.biopython_p19_phylo_functions.run_paml()`

- **作用：** 对序列文件运行 PAML 程序 `codeml`、`baseml` 或 `yn00`。
- **处理方法：** 把比对（和树）复制到新的工作目录，用 Biopython 的选项表（`Bio.Phylo.PAML.<program>`）写出控制文件，以控制文件为唯一参数运行程序并解析结果文件。仅含 Newick 的树会补上 PAML 需要的 `<类群数> 1` 头行。需要 PAML 在 `PATH` 中（`conda install -c bioconda paml`）或指定 `executable`。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.run_paml(program, alignment_path, working_dir, *, tree_path=None, options=None, executable=None, timeout_seconds=3600.0)`。
- **输入：** `program`、PHYLIP/PAML 格式的 `alignment_path`、新的 `working_dir`、可选的 `tree_path`（`yn00` 不用）、按 PAML 名称给出的 `options`（未知名称被拒绝）、`executable` 和 `timeout_seconds`。
- **输出：** `status`、程序 `version`、`return_code`、`results_file`、工作目录中的 `files` 以及解析后的 `results`。
- **示例代码：**

```python
import os, tempfile
from dnakit.phylogeny.biopython_p19_phylo_functions import run_paml

work = os.path.join(tempfile.mkdtemp(), "yn00")
result = run_paml("yn00", "codons.phy", work, options={"icode": 0, "verbose": 0})
print(result["status"], result["version"], result["return_code"])
print(sorted(result["results"]))
print(result["results"]["a"]["b"]["YN00"])
```

- **示例输出：**

```text
completed 4.10.10 0
['a', 'b', 'c', 'd']
{'S': 165.2, 'N': 434.8, 't': 0.0817, 'kappa': 1.1311, 'omega': 0.4805, 'dN': 0.021, 'dN SE': 0.007, 'dS': 0.0437, 'dS SE': 0.0167}
```

## 21) `dnakit.phylogeny.biopython_p19_phylo_functions.read_paml_results()`

- **作用：** 解析 `codeml`、`baseml` 或 `yn00` 的 PAML 结果文件。
- **处理方法：** 调用 `Bio.Phylo.PAML.<program>.read`（`codeml.read`、`baseml.read`、`yn00.read`）并把结果转为普通 JSON 类型。`yn00` 给出每对序列的 Nei-Gojobori、Yang-Nielsen、LWL85 和 LPB93 表，`codeml` 给出位点模型结果，`baseml` 给出模型参数和树长。
- **API：** `dnakit.phylogeny.biopython_p19_phylo_functions.read_paml_results(program, results_file)`。
- **输入：** `program` 和结果文件路径。
- **输出：** 解析后的 `results`。
- **示例代码：**

```python
from dnakit.phylogeny.biopython_p19_phylo_functions import read_paml_results

result = read_paml_results("yn00", "pair_yn00.out")
print(result["results"]["a"]["d"]["YN00"])
print(result["results"]["a"]["d"]["NG86"])
```

- **示例输出：**

```text
{'S': 163.2, 'N': 436.8, 't': 0.2998, 'kappa': 1.1043, 'omega': 1.2262, 'dN': 0.1052, 'dN SE': 0.0164, 'dS': 0.0858, 'dS SE': 0.0242}
{'omega': 1.1566, 'dN': 0.1033, 'dS': 0.0894}
```
