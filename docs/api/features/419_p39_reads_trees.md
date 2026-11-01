# 读段处理与系统树

## 作用

为 P39/P40 整合批次提供明确的 DNA 操作接口，返回可序列化结果。

## 处理方法

函数调用其对应的上游操作；显式指定本地文件，新输出路径不得已存在。UCSC/黑名单下载仅在主动调用时联网。验证环境缺少 NUPACK 与 BLAST+，相关科学计算流程尚未实测。Enformer 使用与架构和适配器配置一致的本地权重。

## API

`dnakit.fasta_fastq.cutadapt_p39_functions`

| Function | Operation |
| --- | --- |
| `cutadapt_filter_reads` | Filter reads by N content, expected errors, average error rate, adapter presence or CASAVA chastity. |
| `cutadapt_rewrite_read_names` | Rewrite read headers using a rename template, length tag or stripped suffix. |
| `cutadapt_clamp_negative_qualities` | Clamp negative Phred qualities to zero with cutadapt --zero-cap. |
| `cutadapt_interleaved_reads` | Read and write interleaved paired-end FASTQ, optionally trimming each mate's 3' adapter. |
| `cutadapt_demultiplex_reads` | Single-index or combinatorial dual-index demultiplexing into a new directory. |

`dnakit.phylogeny.toytree_p39_functions`

| Function | Operation |
| --- | --- |
| `toytree_rf_distance` | Robinson-Foulds bipartition distance between two Newick trees with the same tips. |
| `toytree_consensus` | Majority-rule consensus; min_frequency=0.5 keeps splits present in at least half the trees. |
| `toytree_neighbor_joining` | Infer a neighbor-joining tree from labelled pairwise distances. |
| `toytree_upgma` | Infer a UPGMA tree from labelled pairwise distances. |
| `toytree_mrca` | Return the index and descendant tip labels of the most recent common ancestor. |
| `toytree_monophyly` | Test whether selected tips form a monophyletic group in a rooted Newick tree. |
| `toytree_unique_topologies` | Count distinct rooted or unrooted topologies in a Newick tree sample. |
| `toytree_local_genealogy` | Extract one local ToyTree genealogy from a local tskit file, by tree index or genomic site. |
| `toytree_k80_matrix` | K80 nucleotide transition matrix in toytree's upstream base order. |
| `toytree_tn93_matrix` | TN93 nucleotide transition matrix using upstream parameter and base ordering. |

## 输入

DNA 字符串、JSON 字典/列表、本地基因组文件或模型权重。区间坐标、热力学参数见各函数 docstring；NUPACK 默认材料为 DNA。全部接口支持关键字参数 `verbose`。

## 输出

含完成状态、后端信息和操作结果的字典；文件操作返回输出路径，绘图返回图片路径及对应数据。缺少可选依赖时抛出 `BackendUnavailableError`。各批 manifest 单独记录未验收流程。

**示例代码：**

```python
from dnakit.phylogeny.toytree_p39_functions import toytree_monophyly
print(toytree_monophyly('((a,b),c);', ['a', 'b'], verbose=False)['monophyletic'])
```

**示例输出：**

```text
True
```
