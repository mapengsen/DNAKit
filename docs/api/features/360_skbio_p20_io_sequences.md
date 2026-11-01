# scikit-bio 文件格式、序列查询与比对辅助函数

这一组函数位于 `dnakit.sequence.skbio_p20_io_functions`，是对 [scikit-bio](https://scikit.bio/) 的适配器：QSEQ、Clustal、PHYLIP、Newick、距离矩阵、排序结果、QIIME 2 元数据、BIOM 与 NCBI taxdump 文件的读写，序列查询（屏蔽、位置、基序、反向互补）和比对辅助函数（一致序列、保守性、缺口频率、CIGAR 路径）。scikit-bio 是可选依赖（`dnakit[bio-skibio]`；BIOM 还需 `biom-format`），只有调用函数时才导入。序列必须是 scikit-bio 要求的大写 IUPAC 文本；位置从 0 开始，区间为半开区间。测试版本 scikit-bio 0.7.4。

## 1) `dnakit.sequence.skbio_p20_io_functions.read_skbio_file()`

- **作用：** 用 scikit-bio 读取一个生物信息学文件，并以普通列表和字典返回内容。
- **处理方法：** `file_format` 选择 scikit-bio 读取器：`qseq`（Illumina QSEQ 读段）、`clustal` 与 `phylip`（比对）、`newick`（树）、`lsmat`、`phylip_dm`、`binary_dm`（距离矩阵）、`ordination`、`sample_metadata`（QIIME 2 TSV）、`biom`、`taxdump`（NCBI `*.dmp`，`taxdump_scheme` 指定列名）。一个函数合并了工作清单中的 11 个读取行。列表最多返回 `max_items` 项。
- **API：** `dnakit.sequence.skbio_p20_io_functions.read_skbio_file(path, file_format, *, sequence_type='DNA', qseq_variant='illumina1.8', taxdump_scheme='nodes', max_items=1000)`。
- **输入：** 本地文件路径与格式名；`sequence_type`（DNA 或 RNA）决定比对残基类型，`qseq_variant` 为质量编码。
- **输出：** `format`，以及随格式而定的键，如 `alignment`、`newick`、`ids` 与 `matrix`、`records`、`ordination`，和来源信息。
- **示例代码：**

```python
from pathlib import Path
Path("qseq_reads.txt").write_text("M1\t1\t1\t1\t1\t1\t1\t1\tACG.\tIIII\t1\n")
Path("aln.clustal").write_text("CLUSTAL W (1.83)\n\ns1 ACGT-A\ns2 ACGTTA\n\n")
Path("aln.phy").write_text(" 2 6\ns1        ACGT-A\ns2        ACGTTA\n")
Path("tree.nwk").write_text("((a:1,b:2)90:1,c:3);")
Path("dm.phylip").write_text("2\na 0 1\nb 1 0\n")
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file
reads = read_skbio_file("qseq_reads.txt", "qseq", verbose=False)
print(reads["records"][0]["sequence"], reads["records"][0]["quality"])
for name, fmt in (("aln.clustal", "clustal"), ("aln.phy", "phylip")):
    print(fmt, read_skbio_file(name, fmt, verbose=False)["alignment"])
tree = read_skbio_file("tree.nwk", "newick", verbose=False)
print(tree["newick"], tree["tip_count"])
dm = read_skbio_file("dm.phylip", "phylip_dm", verbose=False)
print(dm["ids"], dm["matrix"])
```

- **示例输出：**

```text
ACG. [40, 40, 40, 40]
clustal {'s1': 'ACGT-A', 's2': 'ACGTTA'}
phylip {'s1': 'ACGT-A', 's2': 'ACGTTA'}
((a:1.0,b:2.0)90:1.0,c:3.0); 3
['a', 'b'] [[0.0, 1.0], [1.0, 0.0]]
```

## 2) `dnakit.sequence.skbio_p20_io_functions.write_skbio_distance_matrix()`

- **作用：** 把对称距离矩阵写成 lsmat 文本或二进制（HDF5）文件。
- **处理方法：** 构建 scikit-bio `DistanceMatrix`，用 `lsmat` 或 `binary_dm` 写出器写文件，`file_format` 选择格式。
- **API：** `dnakit.sequence.skbio_p20_io_functions.write_skbio_distance_matrix(distance_matrix, output_path, *, ids=None, file_format='lsmat')`。
- **输入：** 方阵（列表）、可选 `ids`、输出路径。
- **输出：** 写出的路径、ids 与大小。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_distance_matrix
matrix = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
for fmt, name in (("lsmat", "d.lsmat"), ("binary_dm", "d.h5")):
    out = write_skbio_distance_matrix(matrix, name, ids=["a", "b", "c"], file_format=fmt, verbose=False)
    back = read_skbio_file(name, fmt, verbose=False)
    print(out["format"], out["size"], back["ids"], back["matrix"][2])
```

- **示例输出：**

```text
lsmat 3 ['a', 'b', 'c'] [2.0, 1.5, 0.0]
binary_dm 3 ['a', 'b', 'c'] [2.0, 1.5, 0.0]
```

## 3) `dnakit.sequence.skbio_p20_io_functions.write_skbio_ordination()`

- **作用：** 把排序结果（样本坐标与特征值）写成 scikit-bio 的 ordination 文本格式。
- **处理方法：** 用坐标与特征值构建 `OrdinationResults`，再用 `ordination` 写出器写文件。
- **API：** `dnakit.sequence.skbio_p20_io_functions.write_skbio_ordination(sample_coordinates, eigenvalues, output_path, *, sample_ids=None, proportion_explained=None, method='PCoA', long_method_name='Principal Coordinate Analysis')`。
- **输入：** 样本×轴的坐标、每轴一个特征值、可选样本 id 与解释比例、输出路径。
- **输出：** 写出的路径、样本数与轴数。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_ordination
out = write_skbio_ordination([[0.1, 0.2], [0.3, -0.2], [-0.4, 0.0]], [2.0, 1.0], "pcoa.ord",
                             sample_ids=["a", "b", "c"], proportion_explained=[0.6, 0.4], verbose=False)
print(out["sample_count"], out["axes"])
back = read_skbio_file("pcoa.ord", "ordination", verbose=False)["ordination"]
print(back["samples"]["index"], back["proportion_explained"])
```

- **示例输出：**

```text
3 2
['a', 'b', 'c'] {'0': 0.6, '1': 0.4}
```

## 4) `dnakit.sequence.skbio_p20_io_functions.write_skbio_sample_metadata()`

- **作用：** 写出带类型行的 QIIME 2 样本元数据 TSV。
- **处理方法：** 用各列构建 `SampleMetadata`（数字为数值列，字符串为分类列），调用 `SampleMetadata.save`。
- **API：** `dnakit.sequence.skbio_p20_io_functions.write_skbio_sample_metadata(columns, sample_ids, output_path)`。
- **输入：** 列名到每样本一个值的字典、样本 id 和输出路径。
- **输出：** 写出的路径、样本数与列名。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_sample_metadata
out = write_skbio_sample_metadata({"group": ["x", "y"], "depth": [10, 20]}, ["a", "b"], "meta.tsv", verbose=False)
print(out["columns"])
print(open("meta.tsv").read())
back = read_skbio_file("meta.tsv", "sample_metadata", verbose=False)
print(back["column_types"])
```

- **示例输出：**

```text
['group', 'depth']
id	group	depth
#sk:types	categorical	numeric
a	x	10
b	y	20

{'group': 'categorical', 'depth': 'numeric'}
```

## 5) `dnakit.sequence.skbio_p20_io_functions.write_skbio_biom_table()`

- **作用：** 把计数表写成 BIOM（HDF5）文件。
- **处理方法：** 由样本×特征计数构建 `biom.Table`（特征×样本），用 scikit-bio 的 `biom` 写出器写文件。需要 `biom-format` 包。
- **API：** `dnakit.sequence.skbio_p20_io_functions.write_skbio_biom_table(counts, output_path, *, sample_ids=None, feature_ids=None)`。
- **输入：** 计数（每行一个样本）、可选样本与特征 id、输出路径。
- **输出：** 写出的路径、样本数与特征数。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import read_skbio_file, write_skbio_biom_table
out = write_skbio_biom_table([[1, 2], [3, 4]], "table.biom", sample_ids=["s1", "s2"], feature_ids=["f1", "f2"], verbose=False)
print(out["sample_count"], out["feature_count"])
back = read_skbio_file("table.biom", "biom", verbose=False)
print(back["sample_ids"], back["feature_ids"], back["counts_features_by_samples"])
```

- **示例输出：**

```text
2 2
['s1', 's2'] ['f1', 'f2'] [[1.0, 3.0], [2.0, 4.0]]
```

## 6) `dnakit.sequence.skbio_p20_io_functions.lowercase_masked_sequence()`

- **作用：** 把序列中指定区间改写成小写（软屏蔽）。
- **处理方法：** 由 `[start, end)` 区间生成掩码并调用 `skbio.Sequence.lowercase`。
- **API：** `dnakit.sequence.skbio_p20_io_functions.lowercase_masked_sequence(sequence, mask_intervals, *, sequence_type='DNA')`。
- **输入：** 大写 IUPAC 序列与 0 起点半开区间列表。
- **输出：** 屏蔽后的 `sequence`、被屏蔽位置数和长度。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import lowercase_masked_sequence
r = lowercase_masked_sequence("ACGTACGT", [[1, 3], [6, 8]], verbose=False)
print(r["sequence"], r["masked_positions"])
```

- **示例输出：**

```text
AcgTACgt 4
```

## 7) `dnakit.sequence.skbio_p20_io_functions.find_subsequence()`

- **作用：** 查找子序列的第一次出现位置。
- **处理方法：** 调用 `skbio.Sequence.index`；`start` 与 `end` 限定搜索范围。找不到时返回 `found=False` 而不是报错。
- **API：** `dnakit.sequence.skbio_p20_io_functions.find_subsequence(sequence, subsequence, *, start=None, end=None, sequence_type='DNA')`。
- **输入：** 序列、子序列、可选范围。
- **输出：** `found` 与 0 起点的 `index`。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import find_subsequence
print(find_subsequence("ACGTACGT", "GT", verbose=False)["index"])
print(find_subsequence("ACGTACGT", "GT", start=4, verbose=False)["index"])
print(find_subsequence("ACGT", "TT", verbose=False)["found"])
```

- **示例输出：**

```text
2
6
False
```

## 8) `dnakit.sequence.skbio_p20_io_functions.split_contiguous_segments()`

- **作用：** 按掩码把序列切分成连续片段。
- **处理方法：** 封装 `skbio.Sequence.iter_contiguous`。掩码为区间列表，或字符属于 `included_characters` 的所有位置；`invert` 取补集，`min_length` 丢弃短片段。
- **API：** `dnakit.sequence.skbio_p20_io_functions.split_contiguous_segments(sequence, *, included_intervals=None, included_characters=None, min_length=1, invert=False, sequence_type='DNA')`。
- **输入：** 序列，以及 `included_intervals` 与 `included_characters` 中恰好一个。
- **输出：** 含 `start`、`end`、`sequence` 的 `segments`。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import split_contiguous_segments
r = split_contiguous_segments("ACG-TNA", included_characters="-N", invert=True, verbose=False)
for s in r["segments"]:
    print(s["start"], s["end"], s["sequence"])
```

- **示例输出：**

```text
0 3 ACG
4 5 T
6 7 A
```

## 9) `dnakit.sequence.skbio_p20_io_functions.locate_sequence_positions()`

- **作用：** 定位 IUPAC 序列中的缺口、简并或确定位置。
- **处理方法：** `kind` 选择 `GrammaredSequence.gaps`、`degenerates` 或 `definites`；结果还给出序列是否含确定碱基（`has_definites`）。一个函数覆盖工作清单中的 4 行。
- **API：** `dnakit.sequence.skbio_p20_io_functions.locate_sequence_positions(sequence, kind, *, sequence_type='DNA')`。
- **输入：** 序列与 `kind`。
- **输出：** `positions`、以 `[start, end)` 表示的 `runs`、`count`、`length`、`has_definites`。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import locate_sequence_positions
for kind in ("gaps", "degenerates", "definites"):
    r = locate_sequence_positions("ACG-TNA", kind, verbose=False)
    print(kind, r["positions"], r["runs"], r["has_definites"])
print(locate_sequence_positions("NNN", "definites", verbose=False)["has_definites"])
```

- **示例输出：**

```text
gaps [3] [[3, 4]] True
degenerates [5] [[5, 6]] True
definites [0, 1, 2, 4, 6] [[0, 3], [4, 5], [6, 7]] True
False
```

## 10) `dnakit.sequence.skbio_p20_io_functions.find_nucleotide_motifs()`

- **作用：** 查找核苷酸序列中的嘌呤或嘧啶连续段。
- **处理方法：** 调用 `GrammaredSequence.find_motifs`，类型为 `purine-run` 或 `pyrimidine-run`；`ignore_intervals` 跳过位置，`min_length` 过滤短段。
- **API：** `dnakit.sequence.skbio_p20_io_functions.find_nucleotide_motifs(sequence, motif_type, *, min_length=1, ignore_intervals=None, sequence_type='DNA')`。
- **输入：** 序列、基序类型与可选过滤条件。
- **输出：** 含 `start`、`end`、`sequence` 的 `motifs` 与 `count`。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import find_nucleotide_motifs
r = find_nucleotide_motifs("ACGGGGAGGCGGAG", "purine-run", min_length=2, verbose=False)
for m in r["motifs"]:
    print(m["start"], m["end"], m["sequence"])
print(find_nucleotide_motifs("ACGGGGAGGCGGAG", "pyrimidine-run", min_length=2, verbose=False)["count"])
```

- **示例输出：**

```text
2 9 GGGGAGG
10 14 GGAG
0
```

## 11) `dnakit.sequence.skbio_p20_io_functions.is_reverse_complement_pair()`

- **作用：** 判断两条序列是否互为反向互补。
- **处理方法：** 调用 `NucleotideMixin.is_reverse_complement`，支持 IUPAC 简并码。
- **API：** `dnakit.sequence.skbio_p20_io_functions.is_reverse_complement_pair(sequence_a, sequence_b, *, sequence_type='DNA')`。
- **输入：** 两条同类型序列。
- **输出：** `is_reverse_complement`（真或假）。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import is_reverse_complement_pair
print(is_reverse_complement_pair("AACG", "CGTT", verbose=False)["is_reverse_complement"])
print(is_reverse_complement_pair("AACG", "ACGT", verbose=False)["is_reverse_complement"])
```

- **示例输出：**

```text
True
False
```

## 12) `dnakit.sequence.skbio_p20_io_functions.msa_consensus_sequence()`

- **作用：** 多序列比对的多数一致序列。
- **处理方法：** 构建 `TabularMSA` 并调用 `consensus()`；平局按 scikit-bio 的规则处理。
- **API：** `dnakit.sequence.skbio_p20_io_functions.msa_consensus_sequence(alignment, *, sequence_type='DNA')`。
- **输入：** 名称到比对序列的字典（等长）。
- **输出：** `consensus`、`sequence_count`、`alignment_length`。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import msa_consensus_sequence
aln = {"a": "AC-T", "b": "ACGT", "c": "ACGA"}
r = msa_consensus_sequence(aln, verbose=False)
print(r["consensus"], r["sequence_count"], r["alignment_length"])
```

- **示例输出：**

```text
ACGT 3 4
```

## 13) `dnakit.sequence.skbio_p20_io_functions.msa_conservation()`

- **作用：** 比对每一列的保守性。
- **处理方法：** 按 scikit-bio 默认值调用 `TabularMSA.conservation`（逆 Shannon 不确定度，含缺口的列为 null）。`gap_mode` 与 `degenerate_mode` 可改变处理方式。
- **API：** `dnakit.sequence.skbio_p20_io_functions.msa_conservation(alignment, *, metric='inverse_shannon_uncertainty', degenerate_mode='error', gap_mode='nan', sequence_type='DNA')`。
- **输入：** 比对字典和上述选项。
- **输出：** `conservation`，每列一个值。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import msa_conservation
aln = {"a": "AC-T", "b": "ACGT", "c": "ACGA"}
print(msa_conservation(aln, verbose=False)["conservation"])
print(msa_conservation(aln, gap_mode="ignore", verbose=False)["conservation"])
```

- **示例输出：**

```text
[1.0, 1.0, None, 0.5408520829727552]
[1.0, 1.0, 1.0, 0.5408520829727552]
```

## 14) `dnakit.sequence.skbio_p20_io_functions.msa_gap_frequencies()`

- **作用：** 每个比对列或每条序列的缺口数（或比例）。
- **处理方法：** 调用 `TabularMSA.gap_frequencies`；`axis='sequence'` 每列一个值，`axis='position'` 每条序列一个值；`relative` 返回比例。
- **API：** `dnakit.sequence.skbio_p20_io_functions.msa_gap_frequencies(alignment, *, axis='sequence', relative=False, sequence_type='DNA')`。
- **输入：** 比对字典、`axis`、`relative`。
- **输出：** `frequencies`（按序列时还有 `names`）。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import msa_gap_frequencies
aln = {"a": "AC-T", "b": "ACGT", "c": "ACGA"}
print(msa_gap_frequencies(aln, verbose=False)["frequencies"])
r = msa_gap_frequencies(aln, axis="position", relative=True, verbose=False)
print(r["names"], r["frequencies"])
```

- **示例输出：**

```text
[0, 0, 1, 0]
['a', 'b', 'c'] [0.25, 0.0, 0.0]
```

## 15) `dnakit.sequence.skbio_p20_io_functions.parse_cigar_alignment_path()`

- **作用：** 把 CIGAR 字符串解析成成对比对路径。
- **处理方法：** 调用 `PairAlignPath.from_cigar`，返回各段长度、状态、起始偏移、比对长度以及由路径重新生成的 CIGAR。
- **API：** `dnakit.sequence.skbio_p20_io_functions.parse_cigar_alignment_path(cigar, *, starts=None)`。
- **输入：** CIGAR 字符串与可选 `starts`。
- **输出：** `lengths`、`states`、`starts`、`alignment_length`、`cigar`。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import parse_cigar_alignment_path
r = parse_cigar_alignment_path("3M1D2M", verbose=False)
print(r["lengths"], r["states"], r["alignment_length"], r["cigar"])
```

- **示例输出：**

```text
[3, 1, 2] [[0, 2, 0]] 6 3M1D2M
```

## 16) `dnakit.sequence.skbio_p20_io_functions.apply_alignment_path()`

- **作用：** 根据比对路径，把原始序列还原成带缺口的比对。
- **处理方法：** 由成对 `cigar`，或任意条序列的各段长度 `lengths` 与每段 0/1 缺口标记 `states` 构建 `AlignPath`，再调用带 `gap_char` 的 `AlignPath.to_aligned`。
- **API：** `dnakit.sequence.skbio_p20_io_functions.apply_alignment_path(sequences, *, cigar=None, lengths=None, states=None, starts=None, gap_char='-', sequence_type='DNA')`。
- **输入：** 名称到未比对序列的字典，以及 `cigar` 或 `lengths` 加 `states`、`starts`。
- **输出：** `alignment`（名称到带缺口序列）和 `alignment_length`。
- **示例代码：**

```python
from dnakit.sequence.skbio_p20_io_functions import apply_alignment_path
seqs = {"x": "ACGTAC", "y": "ACGAC"}
print(apply_alignment_path(seqs, cigar="3M1D2M", verbose=False)["alignment"])
r = apply_alignment_path(seqs, lengths=[3, 1, 2], states=[[0, 0], [0, 1], [0, 0]], starts=[0, 0], verbose=False)
print(r["alignment"])
```

- **示例输出：**

```text
{'x': 'ACGTAC', 'y': 'ACG-AC'}
{'x': 'ACGTAC', 'y': 'ACG-AC'}
```
