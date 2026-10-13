# gimmemotifs、tangermeme 与 snapATAC2 的调控区域函数

本页记录 `dnakit.regulation` 中 33 个包装函数：gimmemotifs 12 个（PFM 构建与比较、基序扫描与富集统计、聚类、峰与区间处理、启动子窗口）、tangermeme 13 个（模型突变、变异效应、消融、基序间距与设计、seqlet、k-mer、位点切取、背景位点、信息含量与基序放置）和 snapATAC2 8 个（TSS 富集、FRiP、片段长度与染色体计数、基因活性矩阵、区域-基因连接、基序富集与峰合并）。所有函数只读取本地文件、显式序列或调用方提供的 PyTorch 模型，不下载基因组或基序数据库，也不训练模型。依赖为可选包：[gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs)、[snapATAC2](https://github.com/scverse/snapatac2) 与 [tangermeme](https://github.com/jmschrei/tangermeme)（tangermeme 另需 PyTorch）。

依赖缺失时，函数会报错并给出 `bio-gimmemotifs`、`bio-tangermeme` 或 `bio-snapatac2` 的安装提示。上游包只在函数调用时导入。

示例中的导入路径是各模块文件（`dnakit.regulation.<模块>`）；统一导出由后续接入时完成。

## 示例数据

下面的示例共用这些小文件。它们只用于演示接口，计数和得分不代表真实数据。示例中的基因组 FASTA 由示例代码生成。

**示例文件：** `fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	120	330	cellA	2
chr1	600	800	cellA	1
chr2	50	250	cellA	1
chr1	650	900	cellB	1
chr2	50	250	cellB	1
chr2	1200	1400	cellB	3
chr1	2000	2200	cellB	1
```

**示例文件：** `genes.gtf`

```text
chr1	test	gene	90	400	.	+	.	gene_id "G1"; gene_name "GENEA";
chr1	test	transcript	90	400	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENEA";
chr1	test	gene	1000	1300	.	-	.	gene_id "G2"; gene_name "GENEB";
chr1	test	transcript	1000	1300	.	-	.	gene_id "G2"; transcript_id "T2"; gene_name "GENEB";
chr2	test	gene	400	1600	.	+	.	gene_id "G3"; gene_name "GENEC";
chr2	test	transcript	400	1600	.	+	.	gene_id "G3"; transcript_id "T3"; gene_name "GENEC";
```

**示例文件：** `peaks.bed`

```text
chr1	100	400
chr1	1000	1300
chr2	500	900
```

**示例文件：** `motifs.pfm`

```text
>GTTT_demo
0 0 10 0
0 0 0 10
0 0 0 10
0 0 0 10
>ACGT_demo
10 0 0 0
0 10 0 0
0 0 10 0
0 0 0 10
```

**示例文件：** `motifs.meme`

```text
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF GTAT_demo
letter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0
0 0 1 0
0 0 0 1
1 0 0 0
0 0 0 1
```

**示例文件：** `foreground.fa`

```text
>f1
ACGTTGATAAGCCGTAGGATAAGTTCAG
>f2
TTGATAAGGGCCATTAGAAGTCCTAGGA
```

**示例文件：** `background.fa`

```text
>b1
GGGCCCAAATTTGGGCCCAAATTTGGGC
>b2
ATATATCGCGCGTATATAGCGCGTATAT
```

**示例文件：** `peaks.narrowPeak`

```text
chr1	100	400	p1	10	+	200	5.0	3.0	150
chr1	1000	1300	p2	8	+	120	4.0	2.0	100
chr2	500	900	p3	9	+	300	4.5	2.5	200
```

**示例文件：** `summits.bed`

```text
chr1	250	251	summit1	9.5
chr2	600	601	summit2	7.0
```

**示例文件：** `genes6.bed`

```text
chr1	500	900	G1	0	+
chr1	2000	2500	G2	0	-
chr2	400	800	G3	0	+
```

**示例文件：** `loci.bed`

```text
chr1	100	120
chr2	500	520
```

## 1) `dnakit.regulation.motif_from_consensus()`

- **作用：** 由 IUPAC 共识序列生成 PFM 计数矩阵。
- **处理方法：** 调用 gimmemotifs `motif_from_consensus()`；每个位置的总计数为 `n_sites`，简并碱基把计数平均分给对应碱基。
- **API：** `dnakit.regulation.motif_from_consensus(consensus, *, n_sites=1200, verbose=True)`。
- **输入：** IUPAC 核苷酸共识字符串，长度不超过 100。
- **输出：** `consensus`、`length` 和 `pfm`（每行 A、C、G、T 计数）。
- **示例代码：**

```python
from dnakit.regulation import motif_from_consensus

result = motif_from_consensus("GATAAG", n_sites=100)
print(result["consensus"], result["length"])
print(result["pfm"][0], result["pfm"][2])
```

- **示例输出：**

```text
GATAAG 6
[0, 0, 100, 0] [0, 0, 0, 100]
```

## 2) `dnakit.regulation.motif_from_alignment()`

- **作用：** 由等长比对序列生成 PFM 与概率矩阵。
- **处理方法：** 调用 gimmemotifs `motif_from_align()`，`pfm` 为比对中的碱基计数，`ppm` 为 gimmemotifs 给出的概率（含小伪计数）。
- **API：** `dnakit.regulation.motif_from_alignment(sequences, *, verbose=True)`。
- **输入：** 两条以上等长的 DNA 字符串列表。
- **输出：** `consensus`、`length`、`pfm`、`ppm` 和序列数 `sequence_count`。
- **示例代码：**

```python
from dnakit.regulation import motif_from_alignment

result = motif_from_alignment(["ACGT", "ACGA", "ACGT"])
print(result["consensus"], result["length"])
print(result["pfm"][3])
```

- **示例输出：**

```text
ACGw 4
[1, 0, 0, 2]
```

## 3) `dnakit.regulation.compare_motif_positions()`

- **作用：** 逐位置比较两个等长 PFM 的相似性。
- **处理方法：** 两矩阵先按行归一化并加入 gimmemotifs 的伪计数，再逐位置调用 `akl`（平均 KL 相似性）、`chisq` 或 `ssd`（平方距离和）。
- **API：** `dnakit.regulation.compare_motif_positions(pfm_a, pfm_b, *, metric="akl", verbose=True)`。
- **输入：** 两个等长 PFM（每行四个非负计数）和度量名 `metric`。
- **输出：** `position_scores`（逐位置得分）、`mean_score`、`metric` 和 `length`。
- **示例代码：**

```python
from dnakit.regulation import compare_motif_positions

pfm_a = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0]]
pfm_b = [[10, 0, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]]
akl = compare_motif_positions(pfm_a, pfm_b)
ssd = compare_motif_positions(pfm_a, pfm_b, metric="ssd")
print([round(v, 3) for v in akl["position_scores"]])
print([round(v, 3) for v in ssd["position_scores"]])
```

- **示例输出：**

```text
[10.0, -3.815, -3.815]
[2.0, 0.0, 0.0]
```

## 4) `dnakit.regulation.motif_score_correlation()`

- **作用：** 计算两个 PFM 得分之间的最佳 Pearson 相关。
- **处理方法：** 调用 gimmemotifs `seqcor()`，在一个 7-mer de Bruijn 序列上扫描两个基序并比较得分曲线。
- **API：** `dnakit.regulation.motif_score_correlation(pfm_a, pfm_b, *, verbose=True)`。
- **输入：** 两个 PFM。
- **输出：** `score`（相关系数）、`position`（第二个基序的相对偏移）和 `strand`（+1 或 -1）。
- **示例代码：**

```python
from dnakit.regulation import motif_score_correlation

pfm_a = [[0, 0, 10, 0], [0, 0, 0, 10], [10, 0, 0, 0], [0, 0, 0, 10]]
pfm_b = [[10, 0, 0, 0], [0, 10, 0, 0], [0, 0, 10, 0], [0, 0, 0, 10]]
result = motif_score_correlation(pfm_a, pfm_b)
print(round(result["score"], 4), result["position"], result["strand"])
```

- **示例输出：**

```text
0.5012 -2 -1
```

## 5) `dnakit.regulation.scan_best_motif_hits()`

- **作用：** 为每个基序和每条 FASTA 序列找出最佳命中。
- **处理方法：** 调用 gimmemotifs `scan_to_best_match()`，比较两条链上的 log-odds 得分（`ln(ppm / 0.25 + 0.01)` 之和）。`start` 为正链坐标，负链命中同样用正链起点表示。
- **API：** `dnakit.regulation.scan_best_motif_hits(sequences, motifs, *, max_sequences=1000, verbose=True)`。
- **输入：** FASTA 文件路径和 PFM 文件路径；序列只含 A、C、G、T 时结果最清楚。
- **输出：** `hits` 列表，每项含 `sequence`、`motif`、`score`、`start`、`end` 和 `strand`；按 FASTA 顺序和基序 ID 排序。
- **示例代码：**

```python
from dnakit.regulation import scan_best_motif_hits

result = scan_best_motif_hits("foreground.fa", "motifs.pfm")
print(result["sequence_count"], result["motif_count"], result["hit_count"])
for hit in result["hits"]:
    print(hit["sequence"], hit["motif"], hit["start"], hit["end"], hit["strand"])
```

- **示例输出：**

```text
2 2 4
f1 ACGT_demo 0 4 +
f1 GTTT_demo 2 6 +
f2 ACGT_demo 17 21 +
f2 GTTT_demo 2 6 +
```

## 6) `dnakit.regulation.motif_enrichment_stats()`

- **作用：** 比较前景与背景 FASTA 中每个基序的富集统计。
- **处理方法：** 调用 gimmemotifs `calc_stats()`，计算 ROC AUC、PR AUC、FPR 处的召回与富集等统计量；不做 z 分数和 GC 校正，因此不需要基因组。
- **API：** `dnakit.regulation.motif_enrichment_stats(foreground, background, motifs, *, verbose=True)`。
- **输入：** 前景 FASTA、背景 FASTA 和 PFM 文件。
- **输出：** `motifs` 列表，每项含 `motif` 和 `stats`（统计量名到数值的映射）；按基序 ID 排序。
- **示例代码：**

```python
from dnakit.regulation import motif_enrichment_stats

result = motif_enrichment_stats("foreground.fa", "background.fa", "motifs.pfm")
for motif in result["motifs"]:
    print(motif["motif"], round(motif["stats"]["roc_auc"], 4))
```

- **示例输出：**

```text
ACGT_demo_ACGT 0.875
GTTT_demo_GTTT 0.25
```

## 7) `dnakit.regulation.cluster_motif_set()`

- **作用：** 按相似性对 PFM 基序进行层次聚类。
- **处理方法：** 调用 gimmemotifs `cluster_motifs()`，按 `metric` 与 `match` 计算两两得分，合并得分超过 `threshold` 的基序；每个簇给出成员和平均代表基序的 ID。
- **API：** `dnakit.regulation.cluster_motif_set(motifs, *, metric="wic", match="total", threshold=0.95, verbose=True)`。
- **输入：** PFM 文件（至少两个基序）；`metric` 为 wic、pcc、ed、distance 或 chisq；`match` 为 total 或 subtotal。
- **输出：** `clusters` 列表，每项含 `representative` 与 `members`；按第一个成员排序。
- **示例代码：**

```python
from dnakit.regulation import cluster_motif_set

result = cluster_motif_set("motifs.pfm", threshold=0.95)
for cluster in result["clusters"]:
    print(cluster["representative"], cluster["members"])
```

- **示例输出：**

```text
ACGT_demo ['ACGT_demo']
GTTT_demo ['GTTT_demo']
```

## 8) `dnakit.regulation.read_peak_summits()`

- **作用：** 读取 MACS2 summits.bed 或 narrowPeak 文件中的峰。
- **处理方法：** 调用 gimmemotifs `read_peak_file_to_df()`。summits 文件的每个峰为 1 bp，narrowPeak 的 `start` 为峰起点加 summit 偏移量。
- **API：** `dnakit.regulation.read_peak_summits(path, *, max_peaks=100000, verbose=True)`。
- **输入：** summits.bed 或 narrowPeak 文件；其他 BED 文件会被拒绝。
- **输出：** `peaks` 列表（`chrom`、`start`、`end`、`name`、`value`、`log_value` 等列）和 `columns`。
- **示例代码：**

```python
from dnakit.regulation import read_peak_summits

result = read_peak_summits("summits.bed")
for peak in result["peaks"]:
    print(peak["chrom"], peak["start"], peak["end"], peak["name"], round(peak["log_value"], 4))
```

- **示例输出：**

```text
chr1 250 251 summit1 2.3514
chr2 600 601 summit2 2.0794
```

## 9) `dnakit.regulation.summit_intervals_from_narrowpeak()`

- **作用：** 按 narrowPeak summit 生成固定宽度区间。
- **处理方法：** 调用 gimmemotifs `narrowpeak_to_bed()`；区间从 `峰起点 + summit - size // 2` 开始，长度为 `size`。
- **API：** `dnakit.regulation.summit_intervals_from_narrowpeak(path, *, size, verbose=True)`。
- **输入：** narrowPeak 文件（第 10 列为 summit 偏移）和区间宽度 `size`。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end` 和 `columns`。
- **示例代码：**

```python
from dnakit.regulation import summit_intervals_from_narrowpeak

result = summit_intervals_from_narrowpeak("peaks.narrowPeak", size=50)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"])
```

- **示例输出：**

```text
chr1 225 275
chr1 1075 1125
chr2 675 725
```

## 10) `dnakit.regulation.equal_size_intervals()`

- **作用：** 把 BED 区间缩放为以中点为中心的固定宽度。
- **处理方法：** 调用 gimmemotifs `write_equalsize_bedfile()`。`keep_columns=False` 只保留前六列，与上游 `peakfile2bedfile()` 的输出形式一致。
- **API：** `dnakit.regulation.equal_size_intervals(path, *, size, keep_columns=True, verbose=True)`。
- **输入：** BED 文件（至少三列）和宽度 `size`。
- **输出：** `intervals`，每项含 `chrom`、`start`、`end` 和 `columns`（其余列）。
- **示例代码：**

```python
from dnakit.regulation import equal_size_intervals

result = equal_size_intervals("peaks.narrowPeak", size=100, keep_columns=False)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["columns"][0])
```

- **示例输出：**

```text
chr1 200 300 p1
chr1 1100 1200 p2
chr2 650 750 p3
```

## 11) `dnakit.regulation.random_genomic_regions()`

- **作用：** 从基因组 FASTA 随机抽取固定宽度区间。
- **处理方法：** 调用 gimmemotifs `create_random_genomic_bedfile()`；指定 `seed` 时结果可复现，结果按染色体和起点排序。
- **API：** `dnakit.regulation.random_genomic_regions(genome, size, n, *, seed=None, verbose=True)`。
- **输入：** 基因组 FASTA、区间宽度 `size`、区间数 `n` 和可选种子 `seed`。
- **输出：** `intervals` 列表，每项含 `chrom`、`start` 和 `end`。
- **示例代码：**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")
from dnakit.regulation import random_genomic_regions

result = random_genomic_regions("genome.fa", 100, 3, seed=1)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"])
```

- **示例输出：**

```text
chr1 484 584
chr2 854 954
chr2 941 1041
```

## 12) `dnakit.regulation.promoter_regions()`

- **作用：** 从基因 BED 文件抽取启动子窗口。
- **处理方法：** 调用 gimmemotifs `create_promoter_bedfile()`：正链基因取起点上游 `size` bp，负链基因取终点下游 `size` bp。gimmemotifs 使用全局随机数，本函数在调用前后保存并恢复随机状态。
- **API：** `dnakit.regulation.promoter_regions(genes, size, n, *, seed=None, verbose=True)`。
- **输入：** 六列基因 BED 文件（染色体、起点、终点、名称、分数、链）、窗口宽度与数量 `n`。
- **输出：** `intervals` 列表，每项含 `chrom`、`start`、`end` 和 `strand`。
- **示例代码：**

```python
from dnakit.regulation import promoter_regions

result = promoter_regions("genes6.bed", 100, 3, seed=1)
for item in result["intervals"]:
    print(item["chrom"], item["start"], item["end"], item["strand"])
```

- **示例输出：**

```text
chr1 400 500 +
chr1 2500 2600 -
chr2 300 400 +
```

## 13) `dnakit.regulation.saturation_mutagenesis()`

- **作用：** 对窗口内每个位置做全部单碱基替换，得到模型归因。
- **处理方法：** 调用 tangermeme `saturation_mutagenesis()`。每个位置的突变后输出减去原序列输出，在四个碱基之间居中；`hypothetical=False` 时只保留观察到的碱基。
- **API：** `dnakit.regulation.saturation_mutagenesis(model, sequences, *, start=0, end=None, target=None, hypothetical=False, batch_size=32, verbose=True)`。
- **输入：** 已设为 eval 模式的 PyTorch 模型和等长 DNA 序列（名称到序列的映射）。
- **输出：** `attributions`，索引为序列、碱基（A、C、G、T）、位置（仅窗口 `[start, end)`）。
- **示例代码：**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import saturation_mutagenesis

result = saturation_mutagenesis(model, seqs, start=0, end=5)
print(result["shape"], result["sequences"])
print([round(v, 4) for v in result["attributions"][0][0]])
```

- **示例输出：**

```text
[2, 4, 5] ['s1', 's2']
[-0.3383, 0.0, 0.0, 0.0, -0.0]
```

## 14) `dnakit.regulation.substitution_effect()`

- **作用：** 计算单碱基替换前后模型输出的变化。
- **处理方法：** 调用 tangermeme `substitution_effect()`；`delta` 为替换后输出减替换前输出。同一序列的多个替换会一起施加。
- **API：** `dnakit.regulation.substitution_effect(model, sequences, substitutions, *, verbose=True)`。
- **输入：** 模型、序列映射，以及 `{"sequence", "position", "base"}` 形式的替换列表。
- **输出：** `y_before`、`y_after`、`delta`（每个序列一行）。
- **示例代码：**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import substitution_effect

result = substitution_effect(model, seqs, [{"sequence": "s1", "position": 3, "base": "G"}])
print([round(row[0], 4) for row in result["delta"]])
```

- **示例输出：**

```text
[-0.9991, 0.0]
```

## 15) `dnakit.regulation.insertion_effect()`

- **作用：** 计算插入碱基前后模型输出的变化。
- **处理方法：** 调用 tangermeme `insertion_effect()`。插入后序列被截回模型长度：默认从右侧截断，`left=True` 时从左侧截断。
- **API：** `dnakit.regulation.insertion_effect(model, sequences, insertions, *, left=False, verbose=True)`。
- **输入：** 模型、序列映射和 `{"sequence", "position", "base"}` 形式的插入列表（位置为原序列中插入点之前的下标）。
- **输出：** `y_before`、`y_after`、`delta` 和 `left`。
- **示例代码：**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import insertion_effect

result = insertion_effect(model, seqs, [{"sequence": "s2", "position": 4, "base": "C"}])
print([round(row[0], 4) for row in result["y_after"]])
```

- **示例输出：**

```text
[0.7906, -4.9901]
```

## 16) `dnakit.regulation.deletion_effect()`

- **作用：** 计算删除碱基前后模型输出的变化。
- **处理方法：** 调用 tangermeme `deletion_effect()`。序列必须比模型输入长，因为删除后要截回模型长度；`left` 选择保留左侧或右侧窗口。
- **API：** `dnakit.regulation.deletion_effect(model, sequences, deletions, *, left=False, verbose=True)`。
- **输入：** 比模型输入长的等长序列、删除列表 `{"sequence", "position"}`。
- **输出：** `y_before`、`y_after`、`delta`。
- **示例代码：**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import deletion_effect

long_seqs = {"s1": seqs["s1"] + "A", "s2": seqs["s2"] + "A"}
result = deletion_effect(model, long_seqs, [{"sequence": "s1", "position": 5}])
print([round(row[0], 4) for row in result["y_after"]])
```

- **示例输出：**

```text
[-1.1958, -8.5005]
```

## 17) `dnakit.regulation.ablate_region()`

- **作用：** 反复打乱一个区域，观察模型输出的变化。
- **处理方法：** 调用 tangermeme `ablate()`；`shuffle` 为 `dinucleotide`（保留二核苷酸计数）或 `shuffle`（随机置换）。固定 `random_state` 时结果可复现。
- **API：** `dnakit.regulation.ablate_region(model, sequences, *, start, end, n=20, shuffle="dinucleotide", random_state=0, verbose=True)`。
- **输入：** 模型、等长序列、区域 `[start, end)`、打乱次数 `n`。
- **输出：** `y_before`（每个序列一行）和 `y_after`（序列、打乱次数、输出三维）。
- **示例代码：**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import ablate_region

result = ablate_region(model, seqs, start=2, end=14, n=2, random_state=5)
print(len(result["y_after"]), len(result["y_after"][0]))
print([round(v[0], 4) for v in result["y_after"][0]])
```

- **示例输出：**

```text
2 2
[1.9214, -1.6687]
```

## 18) `dnakit.regulation.motif_pair_spacing()`

- **作用：** 测量两个以上基序按不同间距放入序列后的模型输出。
- **处理方法：** 调用 tangermeme `space()`。`spacings` 每行给出相邻基序之间的间隔，`start` 为第一个基序的位置（默认居中）。
- **API：** `dnakit.regulation.motif_pair_spacing(model, sequences, motifs, spacings, *, start=None, verbose=True)`。
- **输入：** 基序列表（ACGT 字符串）、间距列表（每行 `len(motifs) - 1` 个非负整数）和模型。
- **输出：** `y_before` 与 `y_after`（间距行、序列、输出）。
- **示例代码：**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import motif_pair_spacing

result = motif_pair_spacing(model, seqs, ["GG", "TT"], [[0], [2]], start=2)
print(result["spacing_count"], len(result["y_after"][0]))
print([round(v[0], 4) for v in result["y_after"][1]])
```

- **示例输出：**

```text
2 2
[-7.6877, -10.1997]
```

## 19) `dnakit.regulation.design_greedy_substitution()`

- **作用：** 向一条序列中贪心地插入基序，使模型输出接近目标值。
- **处理方法：** 调用 tangermeme `greedy_substitution()`。每轮尝试所有位置和候选基序（可含反向互补），保留使均方误差最小的替换，直到达到 `max_iter` 或误差低于 `tol`。
- **API：** `dnakit.regulation.design_greedy_substitution(model, sequence, target, *, motifs=None, reverse_complement=True, max_iter=10, tol=0.001, verbose=True)`。
- **输入：** 一条 DNA 序列、目标输出列表、可选基序列表。
- **输出：** `designed_sequence`、`input_sequence`、`target` 和 `length`。
- **示例代码：**

```python
import torch

class Linear(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("weights", torch.sin(torch.arange(80, dtype=torch.float32)).reshape(4, 20))

    def forward(self, x, *args):
        return (x * self.weights).sum(dim=(1, 2)).unsqueeze(-1)

model = Linear().eval()
seqs = {"s1": "ACGTTGATAAGCCGTAGGAT", "s2": "TTGATAAGGGCCATTAGAAG"}

from dnakit.regulation import design_greedy_substitution

result = design_greedy_substitution(model, seqs["s1"], [5.0], motifs=["GGG"], max_iter=2)
print(result["designed_sequence"], result["length"])
```

- **示例输出：**

```text
ACGTTGATAAGGGGTAGGGG 20
```

## 20) `dnakit.regulation.call_seqlets()`

- **作用：** 从归因矩阵中调用 seqlet（连续高归因区段）。
- **处理方法：** 调用 tangermeme `recursive_seqlets()`，用经验 p 值筛选区段；结果按序列号、起点和终点排序。
- **API：** `dnakit.regulation.call_seqlets(attributions, *, threshold=0.01, min_seqlet_len=4, max_seqlet_len=25, additional_flanks=0, verbose=True)`。
- **输入：** 每条序列一行、每个位置一个数值的二维列表。
- **输出：** `seqlets` 列表，每项含 `example_index`、`start`、`end`、`attribution` 与 `p_value`。
- **示例代码：**

```python
import numpy as np
from dnakit.regulation import call_seqlets

rng = np.random.RandomState(0)
values = rng.randn(20, 300) * 0.5
values[0, 100:112] = 8.0
result = call_seqlets(values.tolist(), threshold=0.01)
print(result["seqlet_count"])
for item in result["seqlets"][:2]:
    print(item["example_index"], item["start"], item["end"])
```

- **示例输出：**

```text
10
0 97 115
4 151 155
```

## 21) `dnakit.regulation.kmer_counts()`

- **作用：** 统计每条序列中的连续 k-mer，或按位置得分加权求和。
- **处理方法：** 调用 tangermeme `kmers()`。给定 `scores` 时，每个 k-mer 出现处加上它所跨越位置的得分之和；否则每次出现计 1。只返回实际出现的 k-mer。
- **API：** `dnakit.regulation.kmer_counts(sequences, k, *, scores=None, verbose=True)`。
- **输入：** 名称到序列的映射、`k`（1 至 8）和可选的位置得分。
- **输出：** `entries` 列表，每项含 `sequence`、`kmer` 和 `value`；按序列和 k-mer 排序。
- **示例代码：**

```python
from dnakit.regulation import kmer_counts

counts = kmer_counts({"s": "ACGT"}, 2)
print([(item["kmer"], item["value"]) for item in counts["entries"]])
weighted = kmer_counts({"s": "ACGT"}, 2, scores={"s": [1.0, 2.0, 3.0, 4.0]})
print([item["value"] for item in weighted["entries"]])
```

- **示例输出：**

```text
[('AC', 1.0), ('CG', 1.0), ('GT', 1.0)]
[3.0, 5.0, 7.0]
```

## 22) `dnakit.regulation.extract_loci_sequences()`

- **作用：** 从基因组 FASTA 按 BED 位点切取固定宽度序列。
- **处理方法：** 调用 tangermeme `extract_loci()`，窗口以每个位点的中点为中心，宽度为 `in_window`。
- **API：** `dnakit.regulation.extract_loci_sequences(loci, genome, *, in_window=200, out_window=1000, max_jitter=0, verbose=True)`。
- **输入：** BED 位点文件（至少三列）和基因组 FASTA。
- **输出：** `loci` 列表，每项含 `chrom`、`start`、`end` 和 `sequence`。
- **示例代码：**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")

from dnakit.regulation import extract_loci_sequences

result = extract_loci_sequences("loci.bed", "genome.fa", in_window=10, out_window=8)
for item in result["loci"]:
    print(item["chrom"], item["start"], item["end"], item["sequence"])
```

- **示例输出：**

```text
chr1 100 120 GATAAGCCGT
chr2 500 520 AAGGGCCATT
```

## 23) `dnakit.regulation.matched_background_loci()`

- **作用：** 从基因组中抽取与输入位点 GC 含量匹配的背景位点。
- **处理方法：** 调用 tangermeme `extract_matching_loci()`：候选窗口按 GC 分箱，去掉 N 比例超过 `max_n_perc` 的窗口，并按输入位点的分布取样；固定 `random_state` 时可复现。
- **API：** `dnakit.regulation.matched_background_loci(loci, genome, *, in_window=200, out_window=1000, max_n_perc=0.1, gc_bin_width=0.02, random_state=0, verbose=True)`。
- **输入：** BED 位点文件、基因组 FASTA 和窗口宽度。
- **输出：** `loci` 列表（`chrom`、`start`、`end`）、`background_count` 和 `input_locus_count`。
- **示例代码：**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")

from dnakit.regulation import matched_background_loci

result = matched_background_loci("loci.bed", "genome.fa", in_window=40, out_window=20)
print(result["input_locus_count"], result["background_count"])
print(all(item["end"] - item["start"] == 40 for item in result["loci"]))
```

- **示例输出：**

```text
2 2
True
```

## 24) `dnakit.regulation.pwm_information_content()`

- **作用：** 计算 PWM 每个位置的信息含量（比特）。
- **处理方法：** 调用 tangermeme `information_content()`，按 `2 - 熵` 计算每个位置的信息含量。
- **API：** `dnakit.regulation.pwm_information_content(pwm, *, verbose=True)`。
- **输入：** 每行四个概率（A、C、G、T）的列表。
- **输出：** `information_content`（逐位置比特值）和 `total`。
- **示例代码：**

```python
from dnakit.regulation import pwm_information_content

result = pwm_information_content([[0.25, 0.25, 0.25, 0.25], [1.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.0, 0.0]])
print(result["information_content"], result["total"])
```

- **示例输出：**

```text
[0.0, 2.0, 1.0] 3.0
```

## 25) `dnakit.regulation.place_motif()`

- **作用：** 把基序写入每条序列，替换或插入碱基。
- **处理方法：** 调用 tangermeme `substitute()` 或 `insert()`。`mode="insert"` 插入后从右侧截回原长度；`start=None` 时基序居中。
- **API：** `dnakit.regulation.place_motif(sequences, motif, *, start=None, mode="substitute", verbose=True)`。
- **输入：** 名称到序列的映射、ACGT 基序和模式 `mode`。
- **输出：** `sequences` 列表，每项含 `name` 与 `sequence`。
- **示例代码：**

```python
from dnakit.regulation import place_motif

seqs = {"s1": "ACGTTGATAAGCCGTAGGAT"}
print(place_motif(seqs, "GGG", start=5)["sequences"][0]["sequence"])
print(place_motif(seqs, "GGG", start=5, mode="insert")["sequences"][0]["sequence"])
```

- **示例输出：**

```text
ACGTTGGGAAGCCGTAGGAT
ACGTTGGGGATAAGCCGTAG
```

## 26) `dnakit.regulation.tss_enrichment()`

- **作用：** 计算每个细胞的转录起始位点富集（TSSe）。
- **处理方法：** 先把 fragment 导入临时 snapATAC2 文件，再调用 `metrics.tsse()`，TSS 来自 GTF 中的转录本；处理后删除临时文件。
- **API：** `dnakit.regulation.tss_enrichment(fragment_file, chrom_sizes, gene_gtf, *, exclude_chroms=("chrM", "M"), min_num_fragments=1, max_fragments=5000000, max_cells=100000, verbose=True)`。
- **输入：** 五列 fragment 文件、染色体长度、GTF 基因注释。
- **输出：** `cells`（每个细胞的 `barcode` 与 `tsse`）、`library_tsse` 与 `frac_overlap_TSS`。
- **示例代码：**

```python
from dnakit.regulation import tss_enrichment

result = tss_enrichment("fragments.tsv", {"chr1": 3000, "chr2": 2000}, "genes.gtf")
for cell in result["cells"]:
    print(cell["barcode"], round(cell["tsse"], 4))
print(round(result["library_tsse"], 4))
```

- **示例输出：**

```text
cellA 0.0
cellB 0.8658
0.8264
```

## 27) `dnakit.regulation.fragment_in_peaks_fraction()`

- **作用：** 计算每个细胞中落在峰内的 fragment 比例（FRiP）。
- **处理方法：** 调用 snapATAC2 `metrics.frip()`。`normalized=True` 时除以细胞的总 fragment 数；`count_as_insertion=True` 时按 Tn5 插入位点计数。
- **API：** `dnakit.regulation.fragment_in_peaks_fraction(fragment_file, chrom_sizes, peak_file, *, normalized=True, count_as_insertion=False, min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_peaks=100000, verbose=True)`。
- **输入：** 五列 fragment 文件、染色体长度和峰 BED 文件（前三列为坐标）。
- **输出：** `cells` 列表，每项含 `barcode` 与 `frip`。
- **示例代码：**

```python
from dnakit.regulation import fragment_in_peaks_fraction

result = fragment_in_peaks_fraction("fragments.tsv", {"chr1": 3000, "chr2": 2000}, "peaks.bed")
for cell in result["cells"]:
    print(cell["barcode"], round(cell["frip"], 4))
```

- **示例输出：**

```text
cellA 0.5
cellB 0.0
```

## 28) `dnakit.regulation.fragment_size_distribution()`

- **作用：** 统计数据集中 fragment 长度的分布。
- **处理方法：** 调用 snapATAC2 `metrics.frag_size_distr()`；`counts[i]` 为长度 `i` bp 的 fragment 数，超过 `max_recorded_size` 的不记录。
- **API：** `dnakit.regulation.fragment_size_distribution(fragment_file, chrom_sizes, *, max_recorded_size=1000, min_num_fragments=1, max_fragments=5000000, max_cells=100000, verbose=True)`。
- **输入：** 五列 fragment 文件和染色体长度。
- **输出：** `counts` 列表、`fragment_count`（记录的 fragment 总数）和 `max_recorded_size`。
- **示例代码：**

```python
from dnakit.regulation import fragment_size_distribution

result = fragment_size_distribution("fragments.tsv", {"chr1": 3000, "chr2": 2000}, max_recorded_size=300)
print(result["fragment_count"], result["counts"][200], result["counts"][210])
```

- **示例输出：**

```text
8 6 1
```

## 29) `dnakit.regulation.chromosome_fragment_counts()`

- **作用：** 统计每个细胞在每条染色体上的 fragment 数。
- **处理方法：** 调用 snapATAC2 `metrics.summary_by_chrom()`；`mode` 为 `count`（fragment 数）、`sum`（读段计数之和）或 `mean`。
- **API：** `dnakit.regulation.chromosome_fragment_counts(fragment_file, chrom_sizes, *, mode="count", min_num_fragments=1, max_fragments=5000000, max_cells=100000, verbose=True)`。
- **输入：** 五列 fragment 文件、染色体长度和统计方式。
- **输出：** `chromosomes`（按名称排序）和 `cells`（每个细胞的 `counts`，顺序与染色体列表一致）。
- **示例代码：**

```python
from dnakit.regulation import chromosome_fragment_counts

result = chromosome_fragment_counts("fragments.tsv", {"chr1": 3000, "chr2": 2000})
print(result["chromosomes"])
for cell in result["cells"]:
    print(cell["barcode"], cell["counts"])
```

- **示例输出：**

```text
['chr1', 'chr2']
cellA [3.0, 1.0]
cellB [2.0, 2.0]
```

## 30) `dnakit.regulation.gene_activity_counts()`

- **作用：** 按基因窗口统计每个细胞的 fragment 数（基因活性矩阵）。
- **处理方法：** 调用 snapATAC2 `pp.make_gene_matrix()`：每个基因的窗口为起点上游 `upstream` bp、基因体（可选）和终点下游 `downstream` bp，基因来自 GTF。
- **API：** `dnakit.regulation.gene_activity_counts(fragment_file, chrom_sizes, gene_gtf, *, upstream=2000, downstream=0, include_gene_body=True, min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_entries=1000000, verbose=True)`。
- **输入：** 五列 fragment 文件、染色体长度、GTF 基因注释和窗口参数。
- **输出：** `genes`、`entries`（`[细胞序号, 基因序号, 计数]`）、`nonzero_count` 与 `total_count`。
- **示例代码：**

```python
from dnakit.regulation import gene_activity_counts

result = gene_activity_counts("fragments.tsv", {"chr1": 3000, "chr2": 2000}, "genes.gtf", upstream=500)
print(result["genes"], result["total_count"])
print(result["entries"][:3])
```

- **示例输出：**

```text
['GENEA', 'GENEB', 'GENEC'] 5
[[0, 0, 2], [0, 2, 1], [1, 2, 2]]
```

## 31) `dnakit.regulation.regulatory_regions_to_genes()`

- **作用：** 把区域与距离窗口内的基因连接起来。
- **处理方法：** 调用 snapATAC2 `tl.init_network_from_annotation()`。区域位于基因起点上游 `upstream` bp 或终点下游 `downstream` bp 范围内时建立连接，距离由 snapATAC2 计算。
- **API：** `dnakit.regulation.regulatory_regions_to_genes(regions, gene_gtf, *, upstream=250000, downstream=250000, id_type="gene_name", coding_gene_only=True, verbose=True)`。
- **输入：** `chr:start-end` 形式的区域列表和 GTF 基因注释。
- **输出：** `edges` 列表，每项含 `region`、`gene` 和 `distance`；按区域和基因排序。
- **示例代码：**

```python
from dnakit.regulation import regulatory_regions_to_genes

result = regulatory_regions_to_genes(["chr1:100-300", "chr1:1000-1200", "chr2:500-700"], "genes.gtf", upstream=250, downstream=250)
for edge in result["edges"]:
    print(edge["region"], edge["gene"], edge["distance"])
```

- **示例输出：**

```text
chr1:100-300 GENEA 11
chr1:1000-1200 GENEB 99
chr2:500-700 GENEC 101
```

## 32) `dnakit.regulation.region_motif_enrichment()`

- **作用：** 检验区域组中的基序富集，对照背景区域。
- **处理方法：** 调用 snapATAC2 `tl.motif_enrichment()`，序列取自基因组 FASTA；`method` 为 `binomial`、`hypergeometric` 或 None（使用 snapATAC2 默认方法）。
- **API：** `dnakit.regulation.region_motif_enrichment(regions, motifs, genome, *, background=None, method=None, verbose=True)`。
- **输入：** 区域组（组名到区域列表的映射）、MEME 基序文件和基因组 FASTA；可选背景区域列表。
- **输出：** `groups` 列表，每组含 `group` 和 `rows`（每个基序的 log2 倍数变化、p 值与校正 p 值）。
- **示例代码：**

```python
open("genome.fa", "w").write(">chr1\n" + "ACGTTGATAAGCCGTAGGAT" * 100 + "\n" + ">chr2\n" + "TTGATAAGGGCCATTAGAAG" * 100 + "\n")

from dnakit.regulation import region_motif_enrichment

result = region_motif_enrichment({"g1": ["chr1:100-300", "chr1:1000-1200"], "g2": ["chr2:500-700"]}, "motifs.meme", "genome.fa")
for group in result["groups"]:
    print(group["group"], group["rows"][0]["id"], round(group["rows"][0]["p-value"], 4))
```

- **示例输出：**

```text
g1 GTAT_demo 0.3333
g2 GTAT_demo 0.3333
```

## 33) `dnakit.regulation.merge_peak_sets()`

- **作用：** 把各组的峰合并为一个不重叠的峰集。
- **处理方法：** 调用 snapATAC2 `tl.merge_peaks()`。每个峰以 summit 为中心，向两侧各延伸 `half_width` bp，重叠部分由 snapATAC2 处理。
- **API：** `dnakit.regulation.merge_peak_sets(peaks_by_group, chrom_sizes, *, half_width=250, verbose=True)`。
- **输入：** 各组的 MACS2 风格峰行（`chrom`、`start`、`end`、`name`、`score`、`strand`、`signal_value`、`p_value`、`q_value`、`summit`）和染色体长度。
- **输出：** `peaks` 列表，每项含 `chrom`、`start`、`end` 与 `groups`（包含该峰的组）。
- **示例代码：**

```python
from dnakit.regulation import merge_peak_sets

peak = {"chrom": "chr1", "start": 100, "end": 400, "name": "p", "score": 10, "strand": ".", "signal_value": 5.0, "p_value": 3.0, "q_value": 2.0, "summit": 200}
result = merge_peak_sets({"A": [peak], "B": [dict(peak, chrom="chr2", start=60, end=300, summit=100)]}, {"chr1": 3000, "chr2": 2000}, half_width=50)
for item in result["peaks"]:
    print(item["chrom"], item["start"], item["end"], item["groups"])
```

- **示例输出：**

```text
chr1 250 351 ['A']
chr2 110 211 ['B']
```

## 验证方式

测试文件是 `others/tests/unit/test_full_G12.py`，共 37 项，每个包装函数至少一项。主要检查如下：

- gimmemotifs：共识与比对生成的 PFM 计数与手工结果一致；序列扫描的正链与负链命中坐标与手工放置的位点一致；`calc_stats()`、`scan_to_best_match()`、`cluster_motifs()`、`seqcor()` 与直接调用上游的结果一致；启动子窗口坐标与手算一致，并确认调用后全局随机状态被恢复。
- tangermeme：替换的输出变化等于两个权重之差；饱和突变的归因等于观察碱基权重减去四个碱基的平均值；插入、删除、消融、基序间距与设计的结果与直接调用上游一致；k-mer 计数与加权和与手算一致；位点序列与基因组切片一致；信息含量与基序放置与手算一致。
- snapATAC2：TSS 富集、FRiP、片段长度与染色体计数与手算或直接调用上游一致；基因活性矩阵与直接调用 `make_gene_matrix()` 一致；区域-基因距离与直接调用一致；基序富集的 p 值与直接调用一致；峰合并的坐标与手算一致。
- 三个包缺失时，报错信息中都包含对应的安装提示（使用 monkeypatch 模拟）。

## 当前范围

已包装的操作见上文 33 个小节。其余函数记录在模块台账中，主要分为以下几类：

- 已有包装覆盖：gimmemotifs 与 snapATAC2 的 MEME、PFM、TRANSFAC、JASPAR 读取（`read_motif_file`）、snapATAC2 的 `import_fragments`、`add_tile_matrix`、`make_peak_matrix`，以及 tangermeme 的 `predict`、`marginalize`、`deep_lift_shap`。
- 与 dnakit 现有函数相同：`rc`、`as_fasta`、one-hot 编码、反向互补、VCF 读取、双核苷酸打乱、GC 含量等。
- 需要外部程序（de novo 工具、gffread）、网络（JASPAR、snapATAC2 数据集与基因组下载）或额外数据文件（GC 索引、同源数据库、基因组背景）的函数未包装。台账中以 `needs-data:` 开头的理由标注需要背景数据的函数。
- 绘图和报告函数不包装。
- 与 DNA 相关但本批未选入的函数（例如 tangermeme 的 `gapped_kmers`、`tfmodisco_seqlets`、`annotate_seqlets`，snapATAC2 的 `add_tf_binding`、`macs3`、`export_fragments`、`export_coverage`）记录为 `other`，理由见台账。
