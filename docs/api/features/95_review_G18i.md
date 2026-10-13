# G18i：gimmemotifs、tangermeme 与 snapATAC2 调控区域函数复核

本页记录 `dnakit.regulation` 中由复核任务 G18i 新增的 36 个包装函数：gimmemotifs 9 个（峰顶合并、峰覆盖度、基序扫描写入 BED、区域打分、基序活性、ROC/KS 统计、PFM 对齐、FASTA 记录计数）、tangermeme 19 个（序列编辑与打乱、随机序列、PWM 共识、坐标换算、FASTA 写出、gapped k-mer、注释计数、模型消融与边缘化、序列设计、PISA 归因与 seqlet）和 snapATAC2 8 个（片段导出、覆盖度轨道、BAM 转片段、接触与数值导入、10x QC 指标、MACS3 峰、TF 结合边）。

所有函数只读取本地文件、显式序列、矩阵或调用方提供的模型，不下载基因组或数据库。依赖为可选包；缺失时报错并给出 `bio-gimmemotifs`、`bio-tangermeme` 或 `bio-snapatac2` 的安装提示。上游包只在调用时导入。

示例中的导入路径是各模块文件（`dnakit.regulation.<模块>`）；统一导出由后续接入时完成。


## 1) `dnakit.regulation.combine_summit_files()`

**示例文件：** `gm_summits1.bed`

```text
chr1	100	101	p1	9.0
chr2	500	501	p2	7.0
```

**示例文件：** `gm_summits2.bed`

```text
chr1	110	111	p3	8.0
chr1	2000	2001	p4	6.0
```

- **作用：** 把多个 MACS2 风格的峰顶文件合并，在窗口内保留最强的峰顶。
- **处理方法：** 调用 gimmemotifs `combine_peaks()`，以本地基因组 FASTA 作为染色体来源；`window` 决定合并的宽度，`scale_value` 为真时先按文件缩放峰值。
- **API：** `dnakit.regulation.combine_summit_files(summit_files, genome_fasta, *, window=200, scale_value=False, verbose=True)`.
- **输入：** 峰顶文件路径列表（1 至 50 个）、FASTA 基因组路径、窗口宽度（1 至 100000）和缩放标志。
- **输出：** `summit_count` 和 `summits`（按染色体和位置排序的 `chrom`、`start`、`end`）。
- **示例代码：**

```python
from pathlib import Path
Path("gm_genome.fa").write_text(
    ">chr1\n" + "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 70 + "\n"
    + ">chr2\n" + "GGGCCCAAATTTGGGCCCAAATTTGATAAGGGCCCAAATTTGGGCCCAAATTTGGGC" * 60 + "\n"
)
from dnakit.regulation import combine_summit_files

result = combine_summit_files(["gm_summits1.bed", "gm_summits2.bed"], "gm_genome.fa", window=200)
print(result["summit_count"], result["window"])
for item in result["summits"]:
    print(item["chrom"], item["start"], item["end"])
```

- **示例输出：**

```text
3 200
chr1 100 101
chr1 2000 2001
chr2 500 501
```

## 2) `dnakit.regulation.peak_coverage_table()`

**示例文件：** `gm_summits1.bed`

```text
chr1	100	101	p1	9.0
chr2	500	501	p2	7.0
```

- **作用：** 读取每个峰顶附近窗口内 bigWig 轨道的平均覆盖度。
- **处理方法：** 调用 gimmemotifs `coverage_table()`；每个峰扩展为以峰顶为中心的 `window` 宽区间，`log_transform` 为真时取 log(1 + x)，`top` 为正数时保留方差最大的峰。
- **API：** `dnakit.regulation.peak_coverage_table(peak_file, coverage_files, *, window=200, log_transform=True, top=-1, verbose=True)`.
- **输入：** 峰文件路径、1 至 50 个本地 bigWig 路径、窗口宽度、log 变换标志和 `top`（-1 表示全部）。
- **输出：** `tracks`（轨道名）、`peaks`（区间标签）和 `values`（每个峰一行、每个轨道一列）。
- **示例代码：**

```python
import pyBigWig
bw = pyBigWig.open("gm_cov.bw", "w")
bw.addHeader([("chr1", 4000), ("chr2", 3000)])
bw.addEntries(["chr1"], [0], ends=[4000], values=[2.0])
bw.close()
from dnakit.regulation import peak_coverage_table

result = peak_coverage_table("gm_summits1.bed", ["gm_cov.bw"], window=200)
print(result["tracks"], result["peaks"])
print([round(v, 4) for v in result["values"][0]], [round(v, 4) for v in result["values"][1]])
```

- **示例输出：**

```text
['gm_cov'] ['chr1:0-200', 'chr2:400-600']
[1.0986] [0.0]
```

## 3) `dnakit.regulation.scan_motifs_to_file()`

**示例文件：** `gm_motifs.pfm`

```text
>GATA
0 0 10 0
0 0 0 10
10 0 0 0
0 10 0 0
```

- **作用：** 把超过阈值的基序命中写入 BED 文件，并返回排序后的命中列表。
- **处理方法：** 调用 gimmemotifs `scan_to_file()`，使用 log-odds 得分，不做 z 标准化和 GC 校正；`cutoff` 为 0 到 1 之间的阈值比例，0 保留全部命中。输出文件必须不存在。
- **API：** `dnakit.regulation.scan_motifs_to_file(sequences_fasta, motifs_pfm, output_bed, *, cutoff=0.0, verbose=True)`.
- **输入：** FASTA 序列文件、PFM 文件、尚不存在的输出 BED 路径和阈值比例。
- **输出：** `hit_count` 和 `hits`（每个命中含 `chrom`、`start`、`end`、`motif`、`score`、`strand`）。
- **示例代码：**

```python
from pathlib import Path
Path("gm_genome.fa").write_text(
    ">chr1\n" + "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 70 + "\n"
    + ">chr2\n" + "GGGCCCAAATTTGGGCCCAAATTTGATAAGGGCCCAAATTTGGGCCCAAATTTGGGC" * 60 + "\n"
)
from pathlib import Path
from dnakit.regulation import scan_motifs_to_file

Path("gm_hits.bed").unlink(missing_ok=True)
result = scan_motifs_to_file("gm_genome.fa", "gm_motifs.pfm", "gm_hits.bed", cutoff=0.0)
print(result["hit_count"] > 0, result["output"])
hit = result["hits"][0]
print(hit["chrom"], hit["start"], hit["end"], hit["motif"], round(hit["score"], 3), hit["strand"])
```

- **示例输出：**

```text
True gm_hits.bed
chr1 13 17 GATA -0.4 +
```

## 4) `dnakit.regulation.scan_region_scores()`

**示例文件：** `gm_regions.txt`

```text
chr1:100-300
chr2:400-600
```

- **作用：** 计算每个区域内每个基序的最佳得分，结果是区域乘基序的表格。
- **处理方法：** 调用 gimmemotifs `scan_regionfile_to_table()`，序列来自 FASTA 基因组；不做 z 标准化或 GC 校正，因此不需要 GC 背景。基因组的每条染色体需要长于 1280 个碱基。
- **API：** `dnakit.regulation.scan_region_scores(regions_file, genome_fasta, motifs_pfm, *, verbose=True)`.
- **输入：** 每行一个 `chrom:start-end` 区域的文件、FASTA 基因组和 PFM 文件。
- **输出：** `regions`、`motifs` 和 `scores`（每个区域一行，列与 `motifs` 对应）。
- **示例代码：**

```python
from pathlib import Path
Path("gm_genome.fa").write_text(
    ">chr1\n" + "ACGTTGATAAGCCGTAGGATAAGTTCAGTTGATAAGGGCCATTAGAAGTCCTAGGA" * 70 + "\n"
    + ">chr2\n" + "GGGCCCAAATTTGGGCCCAAATTTGATAAGGGCCCAAATTTGGGCCCAAATTTGGGC" * 60 + "\n"
)
from dnakit.regulation import scan_region_scores

result = scan_region_scores("gm_regions.txt", "gm_genome.fa", "gm_motifs.pfm")
print(result["motifs"], result["regions"])
print([round(v, 4) for v in result["scores"][0]])
```

- **示例输出：**

```text
['GATA'] ['chr1:100-300', 'chr2:400-600']
[-0.4004]
```

## 5) `dnakit.regulation.motif_activity_scores()`

**示例文件：** `gm_labels.tsv`

```text
label	group
r0	1
r1	1
r2	0
r3	0
r4	1
r5	0
```

**示例文件：** `gm_counts.tsv`

```text
region	GATA	TTGA
r0	3	0
r1	2	1
r2	0	0
r3	0	1
r4	4	0
r5	0	2
```

- **作用：** 用超几何检验比较前景和背景区域中的基序计数，得到每个基序的活性值。
- **处理方法：** 调用 gimmemotifs `moap()` 的 `hypergeom` 方法，输入为预先计算的计数表（`scoring="count"`），不再扫描基因组。
- **API：** `dnakit.regulation.motif_activity_scores(labels_file, motif_counts_file, *, method="hypergeom", verbose=True)`.
- **输入：** 带表头、以区域名为索引的标签表（1 为前景，0 为背景）和同样索引的非负整数计数表。
- **输出：** `motifs`（每行一个基序）、`classes`（标签类别）和 `values`（按基序和类别排列的活性值）。
- **示例代码：**

```python
from dnakit.regulation import motif_activity_scores

result = motif_activity_scores("gm_labels.tsv", "gm_counts.tsv")
print(result["classes"], result["motifs"])
for motif, row in zip(result["motifs"], result["values"]):
    print(motif, [round(v, 3) for v in row])
```

- **示例输出：**

```text
['0', '1'] ['GATA', 'TTGA']
GATA [-0.0, 0.699]
TTGA [0.0, -0.0]
```

## 6) `dnakit.regulation.motif_score_metric()`

- **作用：** 计算前景与背景基序得分的一个 ROC、PR 或富集统计量。
- **处理方法：** 按 `metric` 调用 gimmemotifs `rocmetrics` 中的函数：`roc_auc`、`pr_auc`、`max_fmeasure`、`mncp`、`max_enrichment`、`recall_at_fdr`、`matches_at_fpr`、`phyper_at_fpr`、`fraction_fpr`、`score_at_fpr`、`enr_at_fpr`、`roc_auc_xlim` 和 `roc_values`。非有限结果返回 None。
- **API：** `dnakit.regulation.motif_score_metric(foreground, background, *, metric="roc_auc", fpr=0.01, fdr=0.1, minbg=2, verbose=True)`.
- **输入：** 前景和背景得分列表（有限数值）、统计量名称，以及 `fpr`、`fdr`、`minbg` 参数。
- **输出：** `value`（标量统计量），或 `curve`（`roc_values` 的 `fpr` 与 `tpr`）。
- **示例代码：**

```python
from dnakit.regulation import motif_score_metric

roc = motif_score_metric([3, 4, 5, 6], [1, 2, 3, 4], metric="roc_auc")
print(roc["metric"], round(roc["value"], 4))
curve = motif_score_metric([3, 4, 5, 6], [1, 2, 3, 4], metric="roc_values")
print(len(curve["curve"]["fpr"]), len(curve["curve"]["tpr"]))
```

- **示例输出：**

```text
roc_auc 0.875
5 5
```

## 7) `dnakit.regulation.motif_position_ks()`

- **作用：** 对基序命中位置做 Kolmogorov-Smirnov 检验，判断命中是否偏向序列的某一端。
- **处理方法：** 调用 gimmemotifs `ks_pvalue()` 或 `ks_significance()`；`output` 为 `pvalue` 时返回 p 值，为 `significance` 时返回 -log10(p)。
- **API：** `dnakit.regulation.motif_position_ks(foreground_positions, background_positions=None, *, output="pvalue", verbose=True)`.
- **输入：** 前景命中位置，可选的背景位置，以及输出类型。
- **输出：** `output`、`value`（非有限时为 None）以及两组的样本数。
- **示例代码：**

```python
from dnakit.regulation import motif_position_ks

fg = [10, 20, 30, 40, 50]
bg = [100, 110, 120, 130, 140]
pvalue = motif_position_ks(fg, bg, output="pvalue")
print(pvalue["output"], round(pvalue["value"], 4))
significance = motif_position_ks(fg, bg, output="significance")
print(round(significance["value"], 4))
```

- **示例输出：**

```text
pvalue 0.9616
0.017
```

## 8) `dnakit.regulation.equal_length_motif_pair()`

- **作用：** 把两个 PFM 按位置偏移对齐并补齐为等长。
- **处理方法：** 先把计数按行归一化为概率，再调用 gimmemotifs `make_equal_length()`。`pos` 为 `pfm_b` 相对 `pfm_a` 的偏移，缺失位置补 0.25；`truncate` 可删除一侧或两侧的补齐部分。
- **API：** `dnakit.regulation.equal_length_motif_pair(pfm_a, pfm_b, *, pos=0, truncate=None, verbose=True)`.
- **输入：** 两个 PFM（每行 A、C、G、T 四个计数）、整数偏移（-100 至 100）和截断选项。
- **输出：** `length`、补齐后的 `pfm_a` 与 `pfm_b`（概率）。
- **示例代码：**

```python
from dnakit.regulation import equal_length_motif_pair

result = equal_length_motif_pair([[10, 0, 0, 0], [0, 10, 0, 0]], [[0, 0, 10, 0]], pos=0)
print(result["length"])
print(result["pfm_b"])
```

- **示例输出：**

```text
2
[[0.0, 0.0, 1.0, 0.0], [0.25, 0.25, 0.25, 0.25]]
```

## 9) `dnakit.regulation.fasta_record_count()`

**示例文件：** `gm_records.fa`

```text
>a
ACGT
>b
GGCC
>c
TTAA
```

- **作用：** 统计 FASTA 文件中的序列记录数。
- **处理方法：** 调用 gimmemotifs `number_of_seqs_in_file()`，统计以 `>` 开头的记录。
- **API：** `dnakit.regulation.fasta_record_count(path, *, verbose=True)`.
- **输入：** 一
- **输出：** `
- **示例代码：**

```python
from dnakit.regulation import fasta_record_count

result = fasta_record_count("gm_records.fa")
print(result["record_count"], result["file"])
```

- **示例输出：**

```text
3 gm_records.fa
```

## 10) `dnakit.regulation.shuffle_region()`

- **作用：** 在单条序列的指定区域内打乱碱基，区域外保持不变。
- **处理方法：** 调用 tangermeme `ersatz.shuffle()`，在 one-hot 矩阵上打乱 `[start, end)`，重复 `n` 次，并以 `random_state` 固定随机性。
- **API：** `dnakit.regulation.shuffle_region(sequence, start, end, *, n=1, random_state=None, verbose=True)`.
- **输入：** DNA 字符串、0 起点的半开区间 `[start, end)`、副本数 `n`（1 至 1000）和可选随机种子。
- **输出：** `sequences`（`n` 条打乱后的序列）以及 `start`、`end`。
- **示例代码：**

```python
from dnakit.regulation import shuffle_region

result = shuffle_region("ACGTACGTAA", 2, 6, n=2, random_state=0)
for sequence in result["sequences"]:
    print(sequence[:2], sorted(sequence[2:6]), sequence[6:])
```

- **示例输出：**

```text
AC ['A', 'C', 'G', 'T'] GTAA
AC ['A', 'C', 'G', 'T'] GTAA
```

## 11) `dnakit.regulation.randomize_region()`

- **作用：** 用给定的 A、C、G、T 概率替换区域内的碱基。
- **处理方法：** 调用 tangermeme `ersatz.randomize()`；`probs` 归一化后作用于区域内每个位置，默认均匀分布。
- **API：** `dnakit.regulation.randomize_region(sequence, start, end, *, probs=None, n=1, random_state=None, verbose=True)`.
- **输入：** DNA 字符串、区间、四个非负概率（A、C、G、T）、副本数和可选随机种子。
- **输出：** `sequences`、归一化后的 `probs`，以及区间信息。
- **示例代码：**

```python
from dnakit.regulation import randomize_region

result = randomize_region("ACGTACGTAA", 2, 6, probs=[0, 0, 1, 0], random_state=0)
print(result["sequences"][0])
print(result["probs"])
```

- **示例输出：**

```text
ACGGGGGTAA
[0.0, 0.0, 1.0, 0.0]
```

## 12) `dnakit.regulation.delete_region()`

- **作用：** 删除序列中一段碱基并返回较短的序列。
- **处理方法：** 调用 tangermeme `ersatz.delete()`，删除 `[start, end)` 对应的 one-hot 列。
- **API：** `dnakit.regulation.delete_region(sequence, start, end, *, verbose=True)`.
- **输入：** DNA 字符串和半开区间 `[start, end)`。
- **输出：** `sequence`、`length`、`deleted`。
- **示例代码：**

```python
from dnakit.regulation import delete_region

result = delete_region("ACGTACGTAA", 2, 6)
print(result["sequence"], result["length"], result["deleted"])
```

- **示例输出：**

```text
ACGTAA 6 4
```

## 13) `dnakit.regulation.random_one_hot_sequences()`

- **作用：** 按给定碱基概率生成随机 DNA 序列，用于测试和基准。
- **处理方法：** 调用 tangermeme `utils.random_one_hot()`；未指定 `probs` 时四个碱基等概率。
- **API：** `dnakit.regulation.random_one_hot_sequences(count, length, *, probs=None, random_state=None, verbose=True)`.
- **输入：** 序列数（1 至 1000）、长度（1 至 100000）、可选四个概率和随机种子。
- **输出：** `sequences`，以及 `count`、`length` 和归一化后的 `probs`。
- **示例代码：**

```python
from dnakit.regulation import random_one_hot_sequences

print(random_one_hot_sequences(2, 5, probs=[1, 0, 0, 0])["sequences"])
print(random_one_hot_sequences(2, 6, random_state=7)["sequences"])
```

- **示例输出：**

```text
['AAAAA', 'AAAAA']
['ATCGTG', 'GACCGT']
```

## 14) `dnakit.regulation.place_motif_set()`

- **作用：** 按给定间距把多个基序写入同一条序列。
- **处理方法：** 调用 tangermeme `ersatz.multisubstitute()`；第一个基序从 `start` 开始，后续基序在前一个基序结束后经过 `spacing[i]` 个碱基开始，原位替换，长度不变。
- **API：** `dnakit.regulation.place_motif_set(sequence, motifs, spacing, *, start, verbose=True)`.
- **输入：** DNA 字符串、2 至 10 个基序、长度为基序数减一的间距列表和起始位置。
- **输出：** `sequence`、`start` 与 `spacing`。
- **示例代码：**

```python
from dnakit.regulation import place_motif_set

result = place_motif_set("AAAAAAAAAA", ["GG", "TT"], [2], start=1)
print(result["sequence"])
```

- **示例输出：**

```text
AGGAATTAAA
```

## 15) `dnakit.regulation.local_dinucleotide_shuffle()`

- **作用：** 在局部窗口内做二核苷酸打乱，每个窗口保留二核苷酸计数。
- **处理方法：** 调用 tangermeme `ersatz.local_dinucleotide_shuffle()`；序列按 `bin_size` 分箱，每箱内部打乱，等长序列一起处理。
- **API：** `dnakit.regulation.local_dinucleotide_shuffle(sequences, *, n=1, bin_size=1024, min_bin_size=512, random_state=None, verbose=True)`.
- **输入：** 等长 DNA 字符串列表、每序列的打乱份数、分箱大小和随机种子。
- **输出：** `sequences`（每条输入序列对应 `n` 条打乱结果），以及 `bin_size`、`min_bin_size`。
- **示例代码：**

```python
from dnakit.regulation import local_dinucleotide_shuffle

from collections import Counter
seq = "ACGTTGATAAGCCGTAGGAT" * 75
result = local_dinucleotide_shuffle([seq], n=2, bin_size=512, min_bin_size=256, random_state=0)
print(len(result["sequences"]), len(result["sequences"][0]), len(result["sequences"][0][0]))
print(Counter(result["sequences"][0][0]) == Counter(seq))
```

- **示例输出：**

```text
1 2 1500
True
```

## 16) `dnakit.regulation.pwm_consensus_sequence()`

- **作用：** 把 PWM 转换为共识 DNA 序列。
- **处理方法：** 调用 tangermeme `utils.pwm_consensus()`，每行取最大值所在的碱基。
- **API：** `dnakit.regulation.pwm_consensus_sequence(pwm, *, verbose=True)`.
- **输入：** 每行四个非负数（A、C、G、T）的矩阵，各行不能全为零。
- **输出：** `consensus`、`length`。
- **示例代码：**

```python
from dnakit.regulation import pwm_consensus_sequence

print(pwm_consensus_sequence([[0.7, 0.1, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1], [0.1, 0.1, 0.1, 0.7]])["consensus"])
```

- **示例输出：**

```text
AGT
```

## 17) `dnakit.regulation.example_to_sequence_coords()`

- **作用：** 把样本内坐标换算成基因组坐标。
- **处理方法：** 调用 tangermeme `utils.example_to_fasta_coords()`；基因组坐标等于位点起点加样本内坐标，可按 `window` 居中。
- **API：** `dnakit.regulation.example_to_sequence_coords(examples, loci, *, window=None, one_indexed=False, verbose=True)`.
- **输入：** 样本行（`example_idx`、`start`、`end`）和位点行（`chrom`、`start`、`end`）。
- **输出：** `intervals`，每行含 `chrom`、`start`、`end`。
- **示例代码：**

```python
from dnakit.regulation import example_to_sequence_coords

result = example_to_sequence_coords([{"example_idx": 0, "start": 5, "end": 10}], [{"chrom": "chr1", "start": 100, "end": 200}])
print(result["intervals"])
```

- **示例输出：**

```text
[{'chrom': 'chr1', 'start': 105, 'end': 110}]
```

## 18) `dnakit.regulation.one_hot_to_fasta_file()`

- **作用：** 把 DNA 序列写成 FASTA 文件。
- **处理方法：** 先把序列编码为 one-hot，再调用 tangermeme `io.one_hot_to_fasta()` 写文件；输出文件必须不存在。
- **API：** `dnakit.regulation.one_hot_to_fasta_file(sequences, output, *, headers=None, verbose=True)`.
- **输入：** 等长 DNA 字符串列表、输出路径和可选记录名。
- **输出：** `output`、`record_count`、`headers`。
- **示例代码：**

```python
from dnakit.regulation import one_hot_to_fasta_file

from pathlib import Path
Path("tm_out.fa").unlink(missing_ok=True)
result = one_hot_to_fasta_file(["ACGT", "TTGA"], "tm_out.fa", headers=["a", "b"])
print(result["record_count"], result["headers"])
print(open("tm_out.fa").read().strip())
```

- **示例输出：**

```text
2 ['a', 'b']
> a
ACGT
> b
TTGA
```

## 19) `dnakit.regulation.gapped_kmer_features()`

- **作用：** 统计每条序列中的间隔 k-mer（允许中间有空位）。
- **处理方法：** 调用 tangermeme `kmers.gapped_kmers()`，返回非零计数的 `[序列索引, 特征索引, 计数]`。序列需长于最大窗口。
- **API：** `dnakit.regulation.gapped_kmer_features(sequences, *, min_k=2, max_k=2, max_gap=1, max_len=3, max_gkmers=5, verbose=True)`.
- **输入：** 等长 DNA 字符串列表，以及 k、空位和窗口长度限制。
- **输出：** `counts`（三元组列表）和 `feature_count`。
- **示例代码：**

```python
from dnakit.regulation import gapped_kmer_features

result = gapped_kmer_features(["ACGTTGATAAGCCGTAGGAT", "TTGATAAGGGCCATTAGAAG"])
print(result["feature_count"], result["counts"][:4])
```

- **示例输出：**

```text
125 [[0, 8, 2], [0, 9, 2], [0, 17, 2], [0, 21, 2]]
```

## 20) `dnakit.regulation.ablate_annotations()`

- **作用：** 逐个打乱每个注释区间的碱基，并记录模型在打乱前后的输出。
- **处理方法：** 调用 tangermeme `ablate_annotations()`：每个 `(example_idx, start, end)` 区间在其所在序列内打乱 `n` 次，区间内碱基组成保持不变；`random_state` 固定随机数；模型在 CPU 上运行。
- **API：** `dnakit.regulation.ablate_annotations(model, sequences, annotations, *, n=1, random_state=None, verbose=True)`.
- **输入：** 模型（eval 模式的 torch 模块）、等长 DNA 序列列表、注释行列表（`example_idx`、`start`、`end`，0 起始，end 不含）、`n`（1 至 1000）和可选的 `random_state`。
- **输出：** `before`（每个注释一行）、`after`（每个注释、每次打乱一行）、`annotation_count` 和 `n`。
- **示例代码：**

```python
import torch

class CountC(torch.nn.Module):
    def forward(self, x):
        return x[:, 1, :].sum(dim=-1, keepdim=True)
from dnakit.regulation import ablate_annotations

result = ablate_annotations(CountC().eval(), ["ACGTACGTAA", "TTCCGGAATT"], [{"example_idx": 0, "start": 0, "end": 4}, {"example_idx": 1, "start": 2, "end": 6}], n=2, random_state=0)
print(result["before"])
print(result["after"])
```

- **示例输出：**

```text
[[2.0], [2.0]]
[[[2.0], [2.0]], [[2.0], [2.0]]]
```

## 21) `dnakit.regulation.annotate_seqlets()`

- **作用：** 用 TOMTOM 为每个 seqlet 找出最相似的基序，并给出 p 值。
- **处理方法：** 调用 tangermeme `annotate_seqlets()`，底层是 memelite 的 TOMTOM（纯 Python 与 numba 实现，不需要 MEME 可执行文件）；每个 PWM 的每一行先归一化为和为 1，再取 `n_nearest` 个最显著的命中。
- **API：** `dnakit.regulation.annotate_seqlets(sequences, seqlets, motifs, *, n_nearest=1, verbose=True)`.
- **输入：** 等长 DNA 序列列表；seqlet 行（`example_idx`、`start`、`end`）；基序字典，名称映射到每个位置一行、四列（A、C、G、T）的 PWM；`n_nearest`（1 至基序数）。建议至少 25 个基序。
- **输出：** `matches`：每个 seqlet 一个列表，元素含 `motif`（无显著命中时为 None）和 `p_value`；另有 `seqlet_count` 和 `motif_count`。
- **示例代码：**

```python
from dnakit.regulation import annotate_seqlets

import numpy

database = {"GATA": [[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]}
rng = numpy.random.RandomState(7)
for index in range(30):
    database[f"rand{index}"] = rng.dirichlet(numpy.ones(4), size=6).tolist()
result = annotate_seqlets(["TTCCGGAATTGATACCGGTTAACCGG"], [{"example_idx": 0, "start": 10, "end": 14}], database, n_nearest=2)
best = result["matches"][0][0]
print(result["motif_count"], best["motif"], best["p_value"] < 0.001)
```

- **示例输出：**

```text
31 GATA True
```

## 22) `dnakit.regulation.count_seqlet_annotations()`

- **作用：** 统计每个样本中各注释出现的次数。
- **处理方法：** 调用 tangermeme `count_annotations()`，输入 `[example_idx, annotation_idx]` 行；计数固定为 int64，避免上游默认的 uint8 溢出；`dim` 为 None 时返回矩阵，为 0 或 1 时沿对应轴求和。
- **API：** `dnakit.regulation.count_seqlet_annotations(pairs, *, dim=None, verbose=True)`.
- **输入：** 非负整数的 `[example_idx, annotation_idx]` 行列表（1 至 10000 行）和可选的 `dim`（None、0 或 1）。
- **输出：** `counts`（矩阵或一维列表）和 `dim`。
- **示例代码：**

```python
from dnakit.regulation import count_seqlet_annotations

pairs = [[0, 0], [0, 1], [0, 1], [1, 0]]
print(count_seqlet_annotations(pairs)["counts"])
print(count_seqlet_annotations(pairs, dim=1)["counts"])
```

- **示例输出：**

```text
[[1, 2], [1, 0]]
[3, 1]
```

## 23) `dnakit.regulation.pairwise_annotation_counts()`

- **作用：** 统计同一样本中注释对的出现次数；给出 `max_distance` 时按两个注释之间的间隔分箱。
- **处理方法：** `max_distance` 为 None 时调用 tangermeme `pairwise_annotations()`，每个不同的对在每个样本中计一次，结果为对称矩阵；给出整数 `max_distance` 时调用 `pairwise_annotations_spacing()`，结果为 注释×注释×距离 的三维数组，距离为左区间终点到右区间起点之间的碱基数。
- **API：** `dnakit.regulation.pairwise_annotation_counts(annotations, *, max_distance=None, symmetric=True, verbose=True)`.
- **输入：** `max_distance` 为 None 时输入 `[example_idx, annotation_idx]` 行；否则输入含 `example_idx`、`annotation`、`start`、`end` 的映射行，同一样本内区间不得重叠。`symmetric` 为真或假。
- **输出：** `mode`（`pairs` 或 `spacing`）、`max_distance`、`counts`（嵌套列表）和 `row_count`。
- **示例代码：**

```python
from dnakit.regulation import pairwise_annotation_counts

print(pairwise_annotation_counts([[0, 0], [0, 1], [0, 2], [0, 2]])["counts"])
rows = [{"example_idx": 0, "annotation": 0, "start": 0, "end": 4}, {"example_idx": 0, "annotation": 1, "start": 6, "end": 10}]
spaced = pairwise_annotation_counts(rows, max_distance=5)["counts"]
print(spaced[0][1][2], spaced[1][0][2])
```

- **示例输出：**

```text
[[0, 1, 1], [1, 0, 1], [1, 1, 0]]
1 1
```

## 24) `dnakit.regulation.screen_sequence_design()`

- **作用：** 通过随机筛选，找出模型输出最接近目标值的 DNA 序列。
- **处理方法：** 调用 tangermeme `screen()`：每批生成 `batch_size` 条随机 one-hot 序列并打分，保留平方误差最小的 `n_best` 条；最差保留项在 `tol` 之内或达到 `max_iter` 批时停止，因此一定会结束。`squared_error` 由返回序列的模型输出重新计算。
- **API：** `dnakit.regulation.screen_sequence_design(model, length, target, *, n_best=1, tol=1e-3, max_iter=100, batch_size=32, random_state=None, verbose=True)`.
- **输入：** 模型（eval 模式）、序列长度（1 至 1000）、目标向量（每个模型输出一个数）、`n_best`（1 至 16）、`tol`、`max_iter`（1 至 10000）、`batch_size` 和可选的 `random_state`。`n_best` 大于 1 时，若两个候选的损失完全相等，上游会报错，离散输出的模型请保持为 1。
- **输出：** `sequences`（最优在前）、`predictions`、`squared_error`、`length` 和 `n_best`。
- **示例代码：**

```python
import torch

class CountA(torch.nn.Module):
    def forward(self, x):
        return x[:, 0, :].sum(dim=-1, keepdim=True)

from dnakit.regulation import screen_sequence_design

result = screen_sequence_design(CountA().eval(), 6, [0.0], random_state=0)
print(result["sequences"], result["squared_error"])
```

- **示例输出：**

```text
['TTCTCG'] [0.0]
```

## 25) `dnakit.regulation.greedy_marginalize_design()`

- **作用：** 通过贪心放置基序构建序列，使模型输出接近目标值。
- **处理方法：** 调用 tangermeme `greedy_marginalize()`：第一个基序从 `len // 2 - 1` 处开始，之后的基序在已放置基序 `max_spacing` 个碱基的范围内搜索；只有均方误差下降超过 0.001 时才保留，最多 `max_iter` 个。`reverse_complement` 为真时同时尝试反向互补。未覆盖的位置写作 `N`，两端修剪。
- **API：** `dnakit.regulation.greedy_marginalize_design(model, sequence, target, motifs, *, max_spacing=12, reverse_complement=True, max_iter=10, batch_size=32, verbose=True)`.
- **输入：** 模型（eval 模式）、背景序列（长度至少为 `2 * (max_spacing + 最长基序)`）、目标向量、1 至 20 条基序、`max_spacing`（0 至 100）、`reverse_complement`、`max_iter`（1 至 100）和 `batch_size`。
- **输出：** `sequence`（可能含 `N` 的构建序列）和 `length`。
- **示例代码：**

```python
import torch

class CountT(torch.nn.Module):
    def forward(self, x):
        return x[:, 3, :].sum(dim=-1, keepdim=True)

from dnakit.regulation import greedy_marginalize_design

result = greedy_marginalize_design(CountT().eval(), "A" * 60, [4.0], ["TT"], reverse_complement=False, max_spacing=3)
print(result["sequence"], result["length"])
```

- **示例输出：**

```text
TTNNNTT 7
```

## 26) `dnakit.regulation.marginalize_annotations()`

- **作用：** 把每个注释区间插入背景序列的中心，并记录模型输出的前后变化。
- **处理方法：** 调用 tangermeme `marginalize_annotations()`：区间替换每条背景序列中位置 `len // 2 - 区间长度 // 2` 处的碱基。结果先按注释、再按背景序列索引。
- **API：** `dnakit.regulation.marginalize_annotations(model, sequences, background, annotations, *, verbose=True)`.
- **输入：** 模型（eval 模式）、来源等长序列列表、背景等长序列列表、注释行（`example_idx`、`start`、`end`）。区间长度不得超过背景长度。
- **输出：** `insert_start`、`before` 和 `after`；`before` 与 `after` 的形状为 注释 × 背景序列 × 输出。
- **示例代码：**

```python
import torch

class CountC(torch.nn.Module):
    def forward(self, x):
        return x[:, 1, :].sum(dim=-1, keepdim=True)
from dnakit.regulation import marginalize_annotations

result = marginalize_annotations(CountC().eval(), ["TTCCCCTTTT"], ["AAAAAAAAAA"], [{"example_idx": 0, "start": 2, "end": 6}])
print(result["insert_start"], result["before"], result["after"])
```

- **示例输出：**

```text
[3] [[[0.0]]] [[[4.0]]]
```

## 27) `dnakit.regulation.pisa_attribution()`

- **作用：** 用 PISA 计算模型每个输出对序列各碱基的贡献。
- **处理方法：** 调用 tangermeme `pisa()`，在 CPU 上运行：每条序列与 `n_shuffles` 个二核苷酸打乱序列比较，对梯度乘输入取平均。`hypothetical` 为假时，结果乘以观察到的 one-hot 碱基，所以每个位置只有实际碱基非零。调用前后 NumPy 全局随机状态保持不变。
- **API：** `dnakit.regulation.pisa_attribution(model, sequences, *, n_shuffles=20, hypothetical=False, batch_size=32, random_state=None, verbose=True)`.
- **输入：** 模型（eval 模式）、等长 DNA 序列列表（碱基总数不超过 100000）、`n_shuffles`（1 至 200）、`hypothetical`、`batch_size` 和可选的 `random_state`。
- **输出：** `attributions`（序列 × 输出 × 4 × 长度 的嵌套列表）和 `shape`。
- **示例代码：**

```python
import torch
from dnakit.regulation import pisa_attribution

layer = torch.nn.Linear(4 * 30, 2)
with torch.no_grad():
    layer.weight.copy_(torch.linspace(-1.0, 1.0, 2 * 4 * 30).reshape(2, 4 * 30))
    layer.bias.zero_()
model = torch.nn.Sequential(torch.nn.Flatten(), layer).eval()
result = pisa_attribution(model, ["ACGTTGATAAGCCGTAGGATAAGTTCAGTT"], n_shuffles=3, random_state=0)
print(result["shape"])
```

- **示例输出：**

```text
[1, 2, 4, 30]
```

## 28) `dnakit.regulation.tfmodisco_seqlets()`

- **作用：** 按 TF-MoDISco 的阈值规则，从归因分数中提取 seqlet。
- **处理方法：** 调用 tangermeme `tfmodisco_seqlets()`：对长度为 `window_size` 的窗口求和，由零模型选择阈值，使假发现率接近 `target_fdr`；每个 seqlet 在窗口两侧各保留 `flank` 个碱基。`attribution` 为修剪后区间的求和，结果是确定性的。
- **API：** `dnakit.regulation.tfmodisco_seqlets(attributions, *, window_size=21, flank=10, target_fdr=0.2, min_passing_frac=0.03, max_passing_frac=0.2, weak_threshold_for_counting_sign=0.8, verbose=True)`.
- **输入：** 二维归因矩阵（每行一个样本、每列一个位置，行长至少为 `window_size`）以及窗口、侧翼和阈值参数。
- **输出：** `seqlets`（`example_idx`、`start`、`end`（不含）和 `attribution`）以及 `seqlet_count`。
- **示例代码：**

```python
import numpy
from dnakit.regulation import tfmodisco_seqlets

rng = numpy.random.RandomState(0)
attr = rng.normal(0, 0.1, size=(3, 60))
attr[0, 20:41] += 1.0
attr[1, 25:46] += 1.0
attr[2, 10:31] += 1.0
result = tfmodisco_seqlets(attr.tolist())
for row in result["seqlets"]:
    print(row["example_idx"], row["start"], row["end"])
```

- **示例输出：**

```text
0 10 51
1 15 56
2 0 41
```

## 29) `dnakit.regulation.export_fragment_files()`

**示例文件：** `sa_fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	400	600	cellA	1
chr1	700	900	cellA	1
chr2	50	250	cellB	1
chr1	1000	1200	cellB	1
chr2	1200	1400	cellC	1
chr1	2000	2200	cellC	1
```

- **作用：** 把片段文件按细胞组拆分，每个组写出一个五列的片段文件。
- **处理方法：** 用 snapATAC2 导入片段（每个条形码都保留），再调用 `export_fragments()`；输出文件名为 `{prefix}{组名}{suffix}`。每个输出行是一个保存的片段，读数列写为 1。组名必须以字母或数字开头，已存在的输出文件会被拒绝。
- **API：** `dnakit.regulation.export_fragment_files(fragment_file, chrom_sizes, groups, out_dir, *, prefix='', suffix='.bed', verbose=True)`.
- **输入：** 五列制表符分隔的片段文件、染色体长度字典、覆盖文件中全部条形码的条形码到组名字典、已存在的输出目录，以及可选的前缀和后缀。
- **输出：** `group_count` 和 `groups`（每组的 `group`、`file` 和 `fragment_count`，按组名排序）。
- **示例代码：**

```python
import shutil
from pathlib import Path

from dnakit.regulation import export_fragment_files

shutil.rmtree("sa_export_out", ignore_errors=True)
Path("sa_export_out").mkdir()
result = export_fragment_files(
    "sa_fragments.tsv",
    {"chr1": 3000, "chr2": 2000},
    {"cellA": "G1", "cellB": "G2", "cellC": "G1"},
    "sa_export_out",
    prefix="sa_",
)
print(result["group_count"])
for group in result["groups"]:
    print(group["group"], group["file"], group["fragment_count"])
```

- **示例输出：**

```text
2
G1 sa_G1.bed 5
G2 sa_G2.bed 2
```

## 30) `dnakit.regulation.export_coverage_tracks()`

**示例文件：** `sa_fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	400	600	cellA	1
chr1	700	900	cellA	1
chr2	50	250	cellB	1
chr1	1000	1200	cellB	1
chr2	1200	1400	cellC	1
chr1	2000	2200	cellC	1
```

- **作用：** 为每个细胞组写出一个全基因组覆盖轨道，格式为 bedGraph 或 bigWig。
- **处理方法：** 导入片段后调用 `export_coverage()`，计数策略为片段计数，不做归一化。`bin_size` 决定分箱宽度，`max_fragment_length` 为空时不过滤片段长度。输出文件名为 `{prefix}{组名}.bedgraph` 或 `.bw`。
- **API：** `dnakit.regulation.export_coverage_tracks(fragment_file, chrom_sizes, groups, out_dir, *, bin_size=10, normalization=None, output_format='bedgraph', max_fragment_length=None, prefix='', verbose=True)`.
- **输入：** 片段文件、染色体长度字典、条形码到组名字典、已存在的输出目录，以及分箱宽度（1 至 10000000）、归一化方式（None、RPKM、CPM 或 BPM）、输出格式和片段长度上限。
- **输出：** `track_count`、`output_format`、`bin_size`、`normalization` 和 `tracks`（每组的 `group` 与 `file`）。
- **示例代码：**

```python
import shutil
from pathlib import Path

from dnakit.regulation import export_coverage_tracks

shutil.rmtree("sa_coverage_out", ignore_errors=True)
Path("sa_coverage_out").mkdir()
result = export_coverage_tracks(
    "sa_fragments.tsv",
    {"chr1": 3000, "chr2": 2000},
    {"cellA": "G1", "cellB": "G2", "cellC": "G1"},
    "sa_coverage_out",
    bin_size=100,
)
print(result["track_count"], [track["file"] for track in result["tracks"]])
print(Path("sa_coverage_out/G2.bedgraph").read_text(encoding="utf-8").splitlines()[0])
```

- **示例输出：**

```text
2 ['G1.bedgraph', 'G2.bedgraph']
chr1	1000	1200	1
```

## 31) `dnakit.regulation.fragments_from_bam()`

**示例文件：** `sa_reads.sam`

```text
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:3000
@SQ	SN:chr2	LN:2000
r1	99	chr1	101	60	10M	=	391	300	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r1	147	chr1	391	60	10M	=	101	-300	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r2	99	chr1	601	60	10M	=	791	200	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r2	147	chr1	791	60	10M	=	601	-200	AAAAAAAAAA	IIIIIIIIII	CB:Z:AAAA
r3	99	chr2	51	60	10M	=	241	200	AAAAAAAAAA	IIIIIIIIII	CB:Z:CCCC
r3	147	chr2	241	60	10M	=	51	-200	AAAAAAAAAA	IIIIIIIIII	CB:Z:CCCC
```

- **作用：** 把带条形码标签的 BAM 文件转换为五列片段文件，并返回 QC 比例。
- **处理方法：** 调用 `snapATAC2.pp.make_fragment_file()`，使用 CB 标签作为条形码、最小 MAPQ 30，并按 Tn5 位移（起点 +4、终点 -5）生成片段。示例先用 `pysam.view` 把 SAM 文本转为 BAM。
- **API：** `dnakit.regulation.fragments_from_bam(bam_file, output_fragment_file, *, barcode_tag='CB', is_paired=True, min_mapq=30, mitochondrial_chroms=('chrM', 'M'), verbose=True)`.
- **输入：** 成对末端的 BAM 文件、尚不存在的输出片段文件、两字符的 SAM 标签、MAPQ 下限（0 至 255）和线粒体染色体名称列表。
- **输出：** `fragment_count`、`barcode_count`、`qc`（snapATAC2 的比例）和 `output_file`。
- **示例代码：**

```python
from pathlib import Path

import pysam

from dnakit.regulation import fragments_from_bam

Path("sa_reads.bam").unlink(missing_ok=True)
Path("sa_bam_fragments.tsv").unlink(missing_ok=True)
Path("sa_reads.bam").write_bytes(pysam.view("-b", "sa_reads.sam"))
result = fragments_from_bam("sa_reads.bam", "sa_bam_fragments.tsv")
print(result["fragment_count"], result["barcode_count"])
print(Path("sa_bam_fragments.tsv").read_text(encoding="utf-8"), end="")
```

- **示例输出：**

```text
3 2
chr1	104	395	AAAA	1
chr1	604	795	AAAA	1
chr2	54	245	CCCC	1
```

## 32) `dnakit.regulation.import_contact_pairs()`

**示例文件：** `sa_contacts.tsv`

```text
cellA	chr1	10	chr1	40	1
cellA	chr1	120	chr1	950	2
cellB	chr2	300	chr2	700	1
cellB	chr1	50	chr2	800	1
```

- **作用：** 把单细胞染色质接触对分箱后写入 hdf5 格式的 AnnData 文件。
- **处理方法：** 每行为 `barcode, chrom1, pos1, chrom2, pos2, count`，位置为 0 起始且必须落在染色体长度内。调用 `snapATAC2.pp.import_contacts()`，行为细胞、列为分箱对。
- **API：** `dnakit.regulation.import_contact_pairs(contact_file, chrom_sizes, output_h5ad, *, bin_size=500000, sorted_by_barcode=False, verbose=True)`.
- **输入：** 接触对文件、染色体长度字典、必须以 `.h5ad` 结尾且尚不存在的输出路径、分箱宽度（1 至 10000000）和排序标志。
- **输出：** `cell_count`、`feature_count`（分箱对数）、`contact_pair_count`、`contact_count`（计数之和）、`chromosomes` 和 `output_file`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.regulation import import_contact_pairs

Path("sa_contacts.h5ad").unlink(missing_ok=True)
result = import_contact_pairs(
    "sa_contacts.tsv",
    {"chr1": 1000, "chr2": 1000},
    "sa_contacts.h5ad",
    bin_size=500,
)
print(result["cell_count"], result["feature_count"], result["contact_count"])
```

- **示例输出：**

```text
2 16 5
```

## 33) `dnakit.regulation.import_base_pair_values()`

**示例文件：** `sa_cellA.tsv`

```text
chrom	pos	methyl	unmethyl
chr1	10	3	7
chr1	20	0	4
chr2	5	2	0
```

**示例文件：** `sa_cellB.tsv`

```text
chrom	pos	methyl	unmethyl
chr1	10	1	1
```

- **作用：** 把每个细胞一个的碱基级计数表导入为 hdf5 格式的 AnnData 文件，例如 WGBS 甲基化计数。
- **处理方法：** 输入目录只能包含 `.tsv` 文件，每个文件一个细胞，细胞名为文件名去掉后缀。每个文件的表头必须是 `chrom, pos, methyl, unmethyl`，每个位置只能出现一次。调用 `snapATAC2.pp.import_values()`，并把每个细胞的取值数与文件行数核对。
- **API：** `dnakit.regulation.import_base_pair_values(input_dir, chrom_sizes, output_h5ad, *, verbose=True)`.
- **输入：** 只含每细胞一个 `.tsv` 文件的目录、染色体长度字典，以及必须以 `.h5ad` 结尾且尚不存在的输出路径。
- **输出：** `cell_count`、`total_value_count`、`reference_length`、按细胞排序的 `cells`（每个细胞的 `cell` 与 `value_count`）和 `output_file`。
- **示例代码：**

```python
import shutil
from pathlib import Path

from dnakit.regulation import import_base_pair_values

shutil.rmtree("sa_values_dir", ignore_errors=True)
Path("sa_values_dir").mkdir()
Path("sa_values.h5ad").unlink(missing_ok=True)
for name in ("sa_cellA.tsv", "sa_cellB.tsv"):
    shutil.copy(name, "sa_values_dir")
result = import_base_pair_values(
    "sa_values_dir",
    {"chr1": 100, "chr2": 50},
    "sa_values.h5ad",
)
print(result["cell_count"], result["total_value_count"])
for cell in result["cells"]:
    print(cell["cell"], cell["value_count"])
```

- **示例输出：**

```text
2 4
sa_cellA 3
sa_cellB 1
```

## 34) `dnakit.regulation.tenx_atac_metrics()`

**示例文件：** `sa_make_10x_bam.py`

```text
import random

import pysam


def write_bam(path):
    """30 cells with 300 read pairs inside a peak, and 120 background barcodes with 3 pairs."""

    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 3000}, {"SN": "chr2", "LN": 2000}],
    }
    rng = random.Random(7)
    reads = []
    for index in range(30):
        for k in range(300):
            reads.append((f"C{index:03d}_{k}", f"C{index:03d}", 0, rng.randint(900, 1400), rng.randint(160, 240)))
    for index in range(120):
        for k in range(3):
            tid = rng.choice([0, 1])
            reads.append((f"B{index:03d}_{k}", f"B{index:03d}", tid, rng.randint(0, 1700), rng.randint(160, 240)))
    rng.shuffle(reads)
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for name, barcode, tid, start, tlen in reads:
            for read1 in (True, False):
                record = pysam.AlignedSegment()
                record.query_name = name
                record.query_sequence = "A" * 50
                record.flag = 0x1 | 0x2 | (0x40 if read1 else 0x80) | (0 if read1 else 0x10)
                record.reference_id = tid
                record.reference_start = start if read1 else start + tlen - 50
                record.mapping_quality = 60
                record.cigarstring = "50M"
                record.next_reference_id = tid
                record.next_reference_start = start + tlen - 50 if read1 else start
                record.template_length = tlen if read1 else -tlen
                record.query_qualities = pysam.qualitystring_to_array("I" * 50)
                record.set_tag("CB", barcode, value_type="Z")
                out.write(record)
```

**示例文件：** `sa_genes.gtf`

```text
chr1	t	gene	1000	1500	.	+	.	gene_id "G1"; gene_name "GATA1";
chr1	t	transcript	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
chr1	t	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
```

**示例文件：** `sa_peaks.bed`

```text
chr1	900	1400	p1	10	.
```

- **作用：** 对带条形码的 BAM 运行 10x ATAC 质控流程，返回测序、细胞、文库复杂度、比对和靶向指标。
- **处理方法：** 调用 `snapATAC2.pp.recipe_10x_metrics()`：BAM 转片段、导入、用本地 GTF 计算 TSS 富集，并用给定的 BED 峰计算峰相关指标。示例中的 BAM 由 `sa_make_10x_bam.py` 生成。MACS3 峰calling 不在此运行，峰必须由 BED 文件提供。
- **API：** `dnakit.regulation.tenx_atac_metrics(bam_file, output_fragment_file, output_h5ad, chrom_sizes, gene_gtf, peaks_bed, *, barcode_tag='CB', verbose=True)`.
- **输入：** 带 CB 标签的 BAM、两个尚不存在的输出路径（片段文件与 `.h5ad`）、染色体长度字典、本地 GTF 基因注释和 BED 峰文件。
- **输出：** `metrics`（嵌套字典，键与 snapATAC2 配方一致）、`fragment_file` 和 `h5ad_file`。
- **示例代码：**

```python
from pathlib import Path

from sa_make_10x_bam import write_bam

from dnakit.regulation import tenx_atac_metrics

for name in ("sa_10x.bam", "sa_10x_fragments.tsv", "sa_10x.h5ad"):
    Path(name).unlink(missing_ok=True)
write_bam("sa_10x.bam")
result = tenx_atac_metrics(
    "sa_10x.bam",
    "sa_10x_fragments.tsv",
    "sa_10x.h5ad",
    {"chr1": 3000, "chr2": 2000},
    "sa_genes.gtf",
    "sa_peaks.bed",
)
metrics = result["metrics"]
print(metrics["Sequencing"]["Sequenced_read_pairs"], metrics["Targeting"]["Number_of_peaks"], metrics["Cells"]["Number_of_cells"])
```

- **示例输出：**

```text
9360.0 1 31
```

## 35) `dnakit.regulation.call_peaks_macs3()`

- **作用：** 按细胞组用 MACS3 调用开放染色质峰。
- **处理方法：** 把片段导入临时的 snapATAC2 文件，按组调用 `snapATAC2.tl.macs3()`，并把每组的峰按染色体和位置排序。示例中的片段由代码生成，每组峰的第一个峰会被打印。
- **API：** `dnakit.regulation.call_peaks_macs3(fragment_file, chrom_sizes, groups, *, qvalue=0.05, shift=-100, extsize=200, verbose=True)`.
- **输入：** 片段文件、染色体长度字典、覆盖文件中全部条形码的组映射，以及 q 值阈值（0 到 1，不含 0）、位移和扩展长度。
- **输出：** `group_count` 和 `groups`（每组的 `group`、`peak_count` 与 `peaks`，峰含 `chrom`、`start`、`end`、`q_value` 等 MACS3 列）。
- **示例代码：**

```python
import random

from dnakit.regulation import call_peaks_macs3

rng = random.Random(3)
rows = []
for _ in range(60):
    start = 1000 + rng.randint(0, 200)
    rows.append(f"chr1\t{start}\t{start + 150}\tcellA\t1")
for _ in range(60):
    start = 500 + rng.randint(0, 100)
    rows.append(f"chr2\t{start}\t{start + 150}\tcellB\t1")
for _ in range(20):
    start = rng.randint(0, 1800)
    rows.append(f"chr1\t{start}\t{start + 150}\tcellC\t1")
with open("sa_clustered.tsv", "w", encoding="utf-8") as handle:
    handle.write("\n".join(rows) + "\n")
result = call_peaks_macs3(
    "sa_clustered.tsv",
    {"chr1": 3000, "chr2": 2000},
    {"cellA": "G1", "cellB": "G2", "cellC": "G1"},
)
for group in result["groups"]:
    peak = group["peaks"][0]
    print(group["group"], group["peak_count"], peak["chrom"], peak["start"], peak["end"])
```

- **示例输出：**

```text
G1 1 chr1 934 1420
G2 1 chr2 415 834
```

## 36) `dnakit.regulation.add_tf_binding_edges()`

**示例文件：** `sa_network.gtf`

```text
chr1	t	gene	900	1500	.	+	.	gene_id "G1"; gene_name "GATA1";
chr1	t	transcript	900	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
chr1	t	exon	900	1500	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GATA1";
chr2	t	gene	400	800	.	+	.	gene_id "G2"; gene_name "OTHER";
chr2	t	transcript	400	800	.	+	.	gene_id "G2"; transcript_id "T2"; gene_name "OTHER";
chr2	t	exon	400	800	.	+	.	gene_id "G2"; transcript_id "T2"; gene_name "OTHER";
```

**示例文件：** `sa_motifs.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF GATA1
letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0
 0 0 1 0
 1 0 0 0
 0 0 0 1
 1 0 0 0
 1 0 0 0
 0 0 1 0
 0 0 1 0
 1 0 0 0
 0 0 0 1
 0 1 0 0
```

- **作用：** 由 GTF 建立区域到基因的网络，并通过扫描本地基因组为 motif 添加区域结合边。
- **处理方法：** 先用 `snapATAC2.tl.init_network_from_annotation()` 以基因名建立网络（上下游各 250 kb），再调用 `snapATAC2.tl.add_tf_binding()` 在区域序列中以 p 值阈值扫描 MEME motif。返回的是 motif 到区域的边。示例中的基因组由代码生成，motif 位于 chr1:1100-1110。
- **API：** `dnakit.regulation.add_tf_binding_edges(regions, gtf_file, motif_file, genome_fasta, *, pvalue=1e-5, upstream=250000, downstream=250000, coding_gene_only=True, verbose=True)`.
- **输入：** `chrom:start-end` 区域列表、本地 GTF、MEME motif 文件、本地 FASTA 基因组，以及 p 值阈值（0 到 1，不含 0）和上下游距离。
- **输出：** `motif_count`、`region_count`、`gene_count`、`binding_edge_count` 和按 motif 与区域排序的 `binding_edges`。
- **示例代码：**

```python
from pathlib import Path

from dnakit.regulation import add_tf_binding_edges

sequence = "A" * 1100 + "GATAAGGATC" + "A" * 1890
Path("sa_genome.fa").write_text(">chr1\n" + sequence + "\n>chr2\n" + "C" * 2000 + "\n", encoding="utf-8")
result = add_tf_binding_edges(
    ["chr1:1000-1200", "chr1:2000-2200", "chr2:500-700"],
    "sa_network.gtf",
    "sa_motifs.meme",
    "sa_genome.fa",
)
print(result["motif_count"], result["region_count"], result["gene_count"])
for edge in result["binding_edges"]:
    print(edge["motif"], edge["region"])
```

- **示例输出：**

```text
1 3 2
GATA1 chr1:1000-1200
```

## 验证方式

测试文件为 `others/tests/unit/test_review_G18i.py`。每个包装函数至少有一个测试：与上游直接调用的结果一致，或与手工计算一致。缺少依赖时的安装提示也有测试覆盖。随机函数均接受显式随机种子，不修改全局随机状态。


## 当前范围

- gimmemotifs：17 个需要外部数据（genomepy 基因组与注释、GC 索引、同源数据库、OrthoFinder 结果、ROC 报告表或基因组名背景）；2 个需要 MEME 或 HOMER 的 de novo 输出；`run_maelstrom` 的命令行流程；其余为输出或输入辅助函数、别名和绘图数据。`scan_to_file` 的 `fpr` 阈值需要背景文件，因此包装函数只暴露 `cutoff`。
- tangermeme：`chunk`、`unchunk`、`extract_signal` 与 `entropy` 为数组辅助函数，不单独包装；`ersatz.random_one_hot`、`design.screen.random_one_hot`、`space.multisubstitute` 与 `ablate.shuffle` 为重复函数。
- snapATAC2：`link_tf_to_gene`、`prune_network` 和 `pagerank` 只在内存中的调控网络上运行，没有序列或区间步骤，记为 `not-dna`。

