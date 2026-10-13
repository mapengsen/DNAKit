# 可选比对后端：minimap2、edlib、parasail 与 pyhmmer

这一组函数是 `dnakit.alignment` 中基于可选依赖的适配器，分别调用 [mappy](https://github.com/lh3/minimap2)（minimap2 的 Python 绑定）、[edlib](https://github.com/Martinsos/edlib)、[parasail](https://github.com/jeffdaily/parasail) 和 [pyhmmer](https://github.com/althonos/pyhmmer)。它们计算 DNA 序列的映射位置、编辑距离、打分比对和 profile HMM 搜索；不下载参考数据库，也不调用变异检测或组装程序。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

下面的 `reference.fa` 只有两条 DNA 序列，用于演示接口。比对结果只针对这些短序列，不代表真实基因组。

**示例文件：** `reference.fa`

```text
>chr1
CAAGTGTCCGGCCATGTCCCTCTCAGATTCGAACACAGAACAGAGGATGGACACGGGCGGAGCCAGACGAGTCCTTGAAGCGACGCTGAAATTCCGTCATAGATCATGTGAGACTTTCTG
>chr2
GTCGTGTTCAACAGCCGTTGTTCTCTCCAATGGGAGAGGGGCACTTTTTTATCGAGGATACCGCTGTGCAATTAGTGCGTAGCCTAAATCTGGAATACAT
```

## 1) `dnakit.alignment.map_query_sequences_paf()`

- **作用：** 用 minimap2 把一组命名 DNA 查询序列映射到 FASTA 或 `.mmi` 参考上，返回 PAF 风格的记录。
- **处理方法：** 对每条查询调用 `mappy.Aligner.map()`，再把每个命中整理为 PAF 字段。`qstart`、`qend`、`tstart`、`tend` 为 0-based 半开区间；`nmatch` 是匹配碱基数，`alen` 是比对块长度；`edit_distance` 来自 NM 标签，`cigar` 为 minimap2 的 CIGAR。`strand` 为 `-` 时，查询是以反向互补序列比对到参考上的。
- **API：** `dnakit.alignment.map_query_sequences_paf(reference, queries, *, preset="map-ont", max_queries=10000, max_hits_per_query=100)`。
- **输入：** `reference` 是 FASTA 文件或 `build_minimap2_index()` 写出的 `.mmi`；`queries` 是 `{名称: 序列}` 或 DNA 对象，只允许 A、C、G、T、N；`preset` 取 `sr`、`map-ont`、`map-pb`、`map-hifi`、`asm5`、`asm10` 或 `asm20`。查询数超过 `max_queries`，或某条查询的命中数超过 `max_hits_per_query` 时报错。
- **输出：** `records` 列表，每条包含 `qname`、`qlen`、`qstart`、`qend`、`strand`、`tname`、`tlen`、`tstart`、`tend`、`nmatch`、`alen`、`mapq`、`primary`、`edit_distance` 和 `cigar`；另有 `query_count`、`hit_count` 和 `provenance`。
- **示例代码：**

```python
from dnakit.alignment import map_query_sequences_paf

queries = {
    "read1": "TCTCAGATTCGAACACAGAACAGAGGATGGACACGGGCGGAGCCAGACGAGTCCTTGAAG",
    "read2": "ATGACGGAATTTCAGCGTCGCTTCAAGGACTCGTCTGGCTCCGCCCGTGTCCATCCTCTG",
}
result = map_query_sequences_paf("reference.fa", queries)
for record in result["records"]:
    print(record["qname"], record["tname"], record["tstart"], record["tend"], record["strand"], record["cigar"])
```

- **示例输出：**

```text
read1 chr1 20 80 + 60M
read2 chr1 40 100 - 60M
```

## 2) `dnakit.alignment.build_minimap2_index()`

- **作用：** 为 FASTA 参考写出 minimap2 的 `.mmi` 索引，并列出其中的 contig 名称。
- **处理方法：** 调用 `mappy.Aligner(reference, preset=..., fn_idx_out=...)` 写出索引。输出路径已存在时默认报错，只有 `overwrite=True` 才会替换。
- **API：** `dnakit.alignment.build_minimap2_index(reference, output_path, *, preset="map-ont", overwrite=False)`。
- **输入：** FASTA 参考、输出路径（父目录必须已存在）、`preset` 与 `map_query_sequences_paf()` 相同。
- **输出：** `output`（索引路径）、`contig_count`、`contigs`（contig 名称列表）和 `index_bytes`。
- **示例代码：**

```python
from dnakit.alignment import build_minimap2_index

result = build_minimap2_index("reference.fa", "reference.mmi", overwrite=True)
print(result["contig_count"], result["contigs"])
```

- **示例输出：**

```text
2 ['chr1', 'chr2']
```

## 3) `dnakit.alignment.fetch_minimap2_reference_region()`

- **作用：** 取出参考中一个 contig 的指定区间。
- **处理方法：** 通过 mappy 的序列存储切片，坐标为 0-based 半开区间。mappy 对超出 contig 末端的区间会静默截断，所以本函数会检查返回长度，超出时报错。
- **API：** `dnakit.alignment.fetch_minimap2_reference_region(reference, contig, start, end, *, max_reference_bytes=2000000000)`。
- **输入：** FASTA 或 `.mmi` 参考、contig 名称，以及满足 `0 <= start < end` 的整数坐标。参考文件大小超过 `max_reference_bytes` 时报错。
- **输出：** `sequence`、`length`、`contig`、`start` 和 `end`。
- **示例代码：**

```python
from dnakit.alignment import fetch_minimap2_reference_region

result = fetch_minimap2_reference_region("reference.fa", "chr2", 0, 12)
print(result["sequence"], result["length"])
```

- **示例输出：**

```text
GTCGTGTTCAAC 12
```

## 4) `dnakit.alignment.align_dna_edlib()`

- **作用：** 用 Edlib 做全局（`NW`）或嵌入式（`HW`）编辑距离比对，返回扩展 CIGAR 和碱基计数。
- **处理方法：** 调用 `edlib.align(task="path")`。CIGAR 使用 `=`、`X`、`I` 和 `D`：`I` 表示查询中有、目标中没有的碱基，`D` 表示目标中有、查询中没有的碱基。`iupac=True` 时，IUPAC 歧义码与它所代表的每个碱基视为相等；默认情况下，它们只与自身相等。`HW` 模式下，查询整体比对到目标的一个子串，`target_start` 与 `target_end` 给出该子串的位置。
- **API：** `dnakit.alignment.align_dna_edlib(query, target, *, mode="NW", k=-1, iupac=False, max_length=100000)`。
- **输入：** 两个 DNA 序列（允许 IUPAC 字符）；`mode` 为 `NW` 或 `HW`；`k` 为允许的最大编辑距离，`-1` 表示不限制。
- **输出：** `found`、`edit_distance`、`cigar`、`query_start`、`query_end`、`target_start`、`target_end`、`matches`、`mismatches`、`insertions`、`deletions` 和 `identity`。编辑距离超过 `k` 时，`found` 为 `False`，距离和 CIGAR 为 `None`。
- **示例代码：**

```python
from dnakit.alignment import align_dna_edlib

for iupac in (False, True):
    result = align_dna_edlib("ACRT", "ACGT", iupac=iupac)
    print(iupac, result["edit_distance"], result["cigar"])
```

- **示例输出：**

```text
False 1 2=1X1=
True 0 4=
```

## 5) `dnakit.alignment.locate_approximate_pattern()`

- **作用：** 在一段文本中找出编辑距离最低的近似匹配位置，适合短引物或基序的查找。
- **处理方法：** 调用 `edlib.align(mode="HW", task="locations", k=max_edits)`。返回的是达到最低编辑距离的匹配区间，不是所有编辑距离不超过 `max_edits` 的位置。区间为 0-based 半开区间。
- **API：** `dnakit.alignment.locate_approximate_pattern(pattern, text, *, max_edits=0, iupac=False, max_length=100000)`。
- **输入：** DNA 模式序列和文本序列；`max_edits` 为 0 到 1000 之间的整数。
- **输出：** `found`、`edit_distance`、`locations`（每项含 `start` 和 `end`）和 `location_count`。
- **示例代码：**

```python
from dnakit.alignment import locate_approximate_pattern

result = locate_approximate_pattern("ACGT", "TTACGTTTACGTT", max_edits=0)
print(result["edit_distance"], result["locations"])
```

- **示例输出：**

```text
0 [{'start': 2, 'end': 6}, {'start': 8, 'end': 12}]
```

## 6) `dnakit.alignment.calculate_edit_distance_matrix()`

- **作用：** 计算一组命名 DNA 序列两两之间的全局编辑距离矩阵。
- **处理方法：** 对每一对序列调用 `edlib.align(mode="NW", task="distance")`。矩阵是对称的，对角线为 0。
- **API：** `dnakit.alignment.calculate_edit_distance_matrix(sequences, *, iupac=False, max_sequences=100, max_length=5000)`。
- **输入：** `{名称: 序列}` 或 DNA 对象；序列数和每条序列长度都有显式上限。
- **输出：** `names`、`matrix`（按 `names` 顺序排列）和 `pair_count`。
- **示例代码：**

```python
from dnakit.alignment import calculate_edit_distance_matrix

result = calculate_edit_distance_matrix({"a": "ACGT", "b": "ACGT", "c": "ACGA", "d": "TTTT"})
print(result["names"], result["matrix"])
```

- **示例输出：**

```text
['a', 'b', 'c', 'd'] [[0, 0, 1, 3], [0, 0, 1, 3], [1, 1, 0, 4], [3, 3, 4, 0]]
```

## 7) `dnakit.alignment.align_dna_parasail()`

- **作用：** 用 parasail 做全局（`global`）、局部（`local`）或半全局（`semiglobal`）比对，使用替换矩阵和仿射空位罚分。
- **处理方法：** `matrix="custom"` 时，用 `match` 和 `mismatch` 构造 ACGTN 替换矩阵；`matrix="dnafull"` 时使用 parasail 内置的 EDNAFULL 矩阵，此时忽略 `match` 和 `mismatch`。长度为 k 的空位代价是 `gap_open + (k - 1) * gap_extend`。局部和半全局模式的区间只覆盖两条序列都有残基的部分，不含自由的末端空位；`cigar` 则覆盖整个回溯路径，半全局模式下可能包含开头或结尾的自由空位。
- **API：** `dnakit.alignment.align_dna_parasail(query, target, *, mode="global", matrix="custom", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000)`。
- **输入：** 只允许 A、C、G、T 和 N；`mode` 为 `global`、`local` 或 `semiglobal`；`max_cells` 限制 `len(query) * len(target)`。
- **输出：** `score`、`found`、`cigar`、`query_aligned`、`target_aligned`、`query_start`、`query_end`、`target_start`、`target_end`、`matches`、`mismatches`、`insertions`、`deletions` 和 `identity`。
- **示例代码：**

```python
from dnakit.alignment import align_dna_parasail

result = align_dna_parasail(
    "TTTGATTACATTT", "GGGGATTACAGGG", mode="local", match=2, mismatch=-1, gap_open=5, gap_extend=1
)
print(result["score"], result["cigar"], result["query_start"], result["query_end"])
print(result["target_start"], result["target_end"], result["identity"])
```

- **示例输出：**

```text
14 7= 3 10
3 10 1.0
```

## 8) `dnakit.alignment.score_dna_parasail()`

- **作用：** 只计算最优分数和最优格子的终点，不做回溯。
- **处理方法：** 调用 parasail 的 scan 版本（`*_scan_sat`），不保存回溯矩阵，因此比 `align_dna_parasail()` 占用的内存少。局部比对的分数为 0 时，终点返回 `None`。
- **API：** `dnakit.alignment.score_dna_parasail(query, target, *, mode="global", matrix="custom", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_cells=25000000)`。
- **输入：** 与 `align_dna_parasail()` 相同。
- **输出：** `score`、`query_end` 和 `target_end`（均为 0-based 半开坐标）。
- **示例代码：**

```python
from dnakit.alignment import score_dna_parasail

options = {"match": 2, "mismatch": -1, "gap_open": 5, "gap_extend": 1}
result = score_dna_parasail("GATTACA", "TTTGATTACATTT", mode="semiglobal", **options)
print(result["score"], result["query_end"], result["target_end"])
```

- **示例输出：**

```text
14 7 10
```

## 9) `dnakit.alignment.score_dna_parasail_matrix()`

- **作用：** 计算一组命名 DNA 序列两两之间的 parasail 分数矩阵。
- **处理方法：** 只计算上三角，再镜像到下三角。分数与参数由 `score_dna_parasail()` 决定。
- **API：** `dnakit.alignment.score_dna_parasail_matrix(sequences, *, mode="global", matrix="custom", match=5, mismatch=-4, gap_open=10, gap_extend=1, max_sequences=50, max_length=2000, max_cells=25000000)`。
- **输入：** `{名称: 序列}` 或 DNA 对象；序列数、单条序列长度和每对的 DP 格子数都有显式上限。
- **输出：** `names`、`scores`（按 `names` 顺序排列）和 `pair_count`。
- **示例代码：**

```python
from dnakit.alignment import score_dna_parasail_matrix

sequences = {"x": "ACGTACGT", "y": "ACGTTACGT", "z": "TTTTGGGG"}
result = score_dna_parasail_matrix(
    sequences, mode="local", match=2, mismatch=-1, gap_open=5, gap_extend=1
)
print(result["names"], result["scores"])
```

- **示例输出：**

```text
['x', 'y', 'z'] [[16, 11, 2], [11, 18, 4], [2, 4, 16]]
```

## 10) `dnakit.alignment.build_dna_hmm()`

- **作用：** 从一个对齐的 DNA 文件构建 HMMER3 profile HMM，并可选写入文件。
- **处理方法：** 用 pyhmmer 以 DNA 字母表读取多序列比对，再用默认设置的 `Builder` 和 `Background` 构建 HMM。HMM 名称默认取文件名（不含后缀），也可以用 `name` 指定。
- **API：** `dnakit.alignment.build_dna_hmm(alignment, *, name=None, output_path=None, overwrite=False, max_sequences=10000, max_columns=100000)`。
- **输入：** 所有序列长度相同的 DNA 比对文件，缺口用 `-` 表示。`output_path` 为可选的 HMMER3 输出路径。
- **输出：** `name`、`length`（模型长度，即比对列数）、`nseq`、`effective_nseq`、`consensus` 和 `output`。
- **示例文件：** `msa.fa`

```text
>seqA
ACGTACGTTAGCAAGGCTTA
>seqB
ACGTACGATAGCAAGGCTTA
>seqC
ACGTTCGTTAGCAAGG-TTA
>seqD
ACGTACGTTAGCTAGGCTTA
```

- **示例代码：**

```python
from dnakit.alignment import build_dna_hmm

result = build_dna_hmm("msa.fa", output_path="msa.hmm", overwrite=True)
print(result["name"], result["length"], result["nseq"], result["consensus"])
```

- **示例输出：**

```text
msa 20 4 AcgtacgttAgcaAggcttA
```

## 11) `dnakit.alignment.search_dna_hmm()`

- **作用：** 用一个 DNA profile HMM 搜索 DNA 序列文件，同时检查正链和反链。
- **处理方法：** 调用 pyhmmer 的 `nhmmer`，使用长序列管线。`strand` 可取 `both`（默认）、`forward`（只搜正链）或 `reverse`（只搜反链）。`evalue` 是报告阈值。输出坐标为 1-based 闭区间，且 `start <= end`，方向由 `strand` 表示。
- **API：** `dnakit.alignment.search_dna_hmm(profile, sequences, *, evalue=10.0, strand="both", max_hits=1000, max_sequences=100000, max_residues=100000000)`。
- **输入：** 只含一个 DNA profile 的 HMMER3 文件（例如 `build_dna_hmm()` 写出的文件），以及 FASTA 或 FASTQ 格式的序列文件。蛋白 profile 会被拒绝。
- **输出：** `hits` 列表（每个报告的结构域含 `target`、`start`、`end`、`strand`、`hmm_start`、`hmm_end`、`score`、`evalue` 和 `sequence_evalue`）、`hit_count` 和 `sequences_searched`。
- **示例文件：** `targets.fa`

```text
>t_fwd
GGGGACGTACGTTAGCAAGGCTTAGGGG
>t_rev
CCCCTAAGCCTTGCTAACGTACGTCCCC
>t_none
TTTTTTTTTTTTTTTTTTTTTTTTTTTT
>t_mut
CCCCACGTACGATAGCAAGGCTTACCCC
```

- **示例代码：**

```python
from dnakit.alignment import search_dna_hmm

result = search_dna_hmm("msa.hmm", "targets.fa")
for hit in result["hits"]:
    print(hit["target"], hit["start"], hit["end"], hit["strand"])
```

- **示例输出：**

```text
t_fwd 5 24 +
t_rev 5 24 -
t_mut 5 24 +
```

## 12) `dnakit.alignment.read_hmm_profile()`

- **作用：** 读取 HMMER3 文件中每个 profile 的名称、字母表、长度、共识序列和序列计数，不做搜索。
- **处理方法：** 用 `pyhmmer.plan7.HMMFile` 逐个读取 profile，并限制 profile 数量。
- **API：** `dnakit.alignment.read_hmm_profile(path, *, max_profiles=1000)`。
- **输入：** HMMER3 文件路径，可以包含多个 profile。
- **输出：** `profiles` 列表（每项含 `name`、`accession`、`description`、`alphabet`、`length`、`consensus`、`nseq` 和 `effective_nseq`）以及 `profile_count`。
- **示例代码：**

```python
from dnakit.alignment import read_hmm_profile

profile = read_hmm_profile("msa.hmm")["profiles"][0]
print(profile["name"], profile["alphabet"], profile["length"], profile["nseq"])
```

- **示例输出：**

```text
msa DNA 20 4
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-mappy]"           # mappy (minimap2): map_query_sequences_paf 等
python -m pip install "dnakit[bio-edlib]"           # edlib: align_dna_edlib 等
python -m pip install "dnakit[bio-parasail]"        # parasail: align_dna_parasail 等
python -m pip install "dnakit[bio-pyhmmer]"         # pyhmmer: build_dna_hmm 等
```

项目地址：[minimap2](https://github.com/lh3/minimap2)、[edlib](https://github.com/Martinsos/edlib)、[parasail](https://github.com/jeffdaily/parasail)、[pyhmmer](https://github.com/althonos/pyhmmer)。

## 验证方式

测试文件是 `others/tests/unit/test_alignment_backends.py`，共 30 项。主要检查如下：

- mappy：精确正链命中的坐标和 MAPQ 与直接调用 mappy 的结果一致；反向互补命中的链方向和坐标与手工计算一致；一个替换位点的 NM 和匹配碱基数与手工计算一致；索引写出后可以重新读取，区间超出 contig 时报错。
- edlib：`4=1D4=` 的 CIGAR 与编辑距离和直接调用 edlib 的结果一致；IUPAC 选项与 `HW` 的目标区间按手工计算；定位函数的区间与直接调用 edlib 的 `locations` 一致；距离矩阵与逐对调用的结果一致。
- parasail：分数 20、11、10、4 和 14 与手工计算一致；局部与半全局的区间与手工结果一致；分数与直接调用 parasail 的结果一致；分数矩阵与逐对计算一致且对称。
- pyhmmer：模型长度和共识序列与直接调用 pyhmmer 的结果一致，序列数与输入一致；三条序列的命中坐标与手工计算一致（正链、反链和一处替换均为 5 至 24）；`strand` 过滤有效；蛋白 profile 被拒绝。
- 每个模块都有一项测试：缺少可选依赖时报出对应的安装提示。另外检查未知参数、超过显式上限、布尔值作为整数等输入错误。

## 当前范围

这一批覆盖 minimap2 的内存查询映射与索引读写、Edlib 的 NW 和 HW 比对与距离矩阵、parasail 的全局、局部和半全局比对与打分、pyhmmer 的 DNA profile HMM 构建、nhmmer 搜索和元数据读取。以下内容暂未包含：

- minimap2 的剪接比对、`cs` 标签和任意命令行参数；映射结果只给出 CIGAR 和 NM。
- Edlib 的 `SHW`（前缀）模式，以及多个 query 对同一 target 的批量比对。
- parasail 的 IUPAC 打分、蛋白替换矩阵（如 BLOSUM）、striped 和 profile 接口。
- pyhmmer 的蛋白搜索（`hmmsearch`、`hmmscan`）、`jackhmmer` 和 `phmmer`，以及 HMM 的 `cutoffs` 设置。
