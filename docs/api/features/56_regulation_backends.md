# 调控基因组可选后端：基序扫描、fragment 计数与模型归因

这一组函数是 `dnakit.regulation` 中基于可选依赖的适配器，分别调用 [gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs)、[snapATAC2](https://github.com/scverse/snapatac2) 和 [tangermeme](https://github.com/jmschrei/tangermeme)。gimmemotifs 负责读取 PFM 基序、在给定序列中扫描位点以及给 k-mer 打分；snapATAC2 负责由五列 fragment 文件计算细胞质控指标、bin 计数矩阵和峰计数矩阵；tangermeme 调用用户提供的 PyTorch 模型做预测、基序边缘化和 DeepLIFT/SHAP 归因。所有函数只处理显式给出的文件、序列或模型，不下载基因组、背景序列或基序数据库，也不训练模型，也不计算基序富集的显著性。

依赖都是可选的，只有调用对应函数时才会导入。缺少包时会报错并给出安装命令。

## 示例数据

下面的示例都使用这三个小文件。它们只用于演示接口，计数和得分不代表真实数据。

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

**示例文件：** `fragments.tsv`

```text
chr1	100	300	cellA	1
chr1	120	330	cellA	2
chr1	600	800	cellA	1
chrM	10	200	cellA	1
chr1	650	900	cellB	1
chr2	50	250	cellB	1
chr2	1200	1400	cellB	3
```

**示例文件：** `peaks.bed`

```text
chr1	50	400
chr1	590	950
chr2	1100	1500
```

## 1) `dnakit.regulation.read_motif_file()`

- **作用：** 读取 PFM、TRANSFAC 或 JASPAR 格式的基序文件，返回每个基序的 ID、长度、共识序列、信息量、得分范围和 PFM 计数。
- **处理方法：** 调用 gimmemotifs 的 `read_motifs()` 解析文件。返回的 `pfm` 保留文件中的计数，可直接传给 `scan_motif_hits()`。文件大小和基序数量都有显式上限。
- **API：** `dnakit.regulation.read_motif_file(path, *, fmt="pfm", max_bytes=10000000, max_motifs=1000)`。
- **输入：** 文本基序文件路径；`fmt` 取 `"pfm"`、`"transfac"` 或 `"jaspar"`。
- **输出：** 含 `motif_count` 和 `motifs` 的字典。每个基序包含 `id`、`length`、`consensus`、`information_content`、`min_score`、`max_score` 和 `pfm`。
- **示例代码：**

```python
from dnakit.regulation import read_motif_file

result = read_motif_file("motifs.pfm")
for motif in result["motifs"]:
    print(motif["id"], motif["consensus"], motif["length"])
```

- **示例输出：**

```text
GTTT_demo GTTT 4
ACGT_demo ACGT 4
```

## 2) `dnakit.regulation.scan_motif_hits()`

- **作用：** 在给定的 DNA 序列中查找 PFM 基序的位点，同时报告正链和反链的命中。
- **处理方法：** 用 gimmemotifs 的 `Motif.scan_all()` 计算对数优势得分，即每个位置 `ln(ppm / 0.25 + 0.01)` 的求和。阈值是每个基序自身得分范围的比例：`cutoff = min + (max - min) * threshold`。默认同时扫描反向互补链。`max_hits_per_sequence` 限制每个基序在每条序列中的命中数，超过时报错而不是截断。
- **API：** `dnakit.regulation.scan_motif_hits(sequences, motifs, *, threshold=0.9, scan_rc=True, max_hits_per_sequence=100, max_sequences=1000, max_sequence_length=100000)`。
- **输入：** `DNA(...)` 或 `{序列名: 序列}`，序列只能包含 A、C、G、T。`motifs` 是 `{基序 ID: PFM 计数矩阵}`，每行对应一个位置，列依次为 A、C、G、T。
- **输出：** `hit_count` 与 `hits` 列表。每个命中包含 `sequence`、`motif`、`start`（0-based，正向坐标）、`end`（开区间）、`strand` 和 `score`。
- **示例代码：**

```python
from dnakit.regulation import read_motif_file, scan_motif_hits

motifs = {motif["id"]: motif["pfm"] for motif in read_motif_file("motifs.pfm")["motifs"]}
result = scan_motif_hits(
    {"promoter": "GGAAACGGTTTAAC", "control": "CCCCCCCC"},
    motifs,
    threshold=0.9,
)
for hit in result["hits"]:
    print(hit["sequence"], hit["motif"], hit["start"], hit["strand"], round(hit["score"], 3))
```

- **示例输出：**

```text
promoter GTTT_demo 2 - 5.554
promoter GTTT_demo 7 + 5.554
```

## 3) `dnakit.regulation.score_motif_kmers()`

- **作用：** 计算等长 k-mer 与一个 PFM 的对数优势得分。
- **处理方法：** 调用 gimmemotifs 的 `Motif.score_kmer()`。同时返回基序的最大和最小得分，便于把得分换算成阈值。
- **API：** `dnakit.regulation.score_motif_kmers(motif, kmers, *, max_kmers=10000)`。
- **输入：** 一个 PFM 计数矩阵（行为位置，列为 A、C、G、T）；k-mer 字符串列表，长度必须与基序长度相同，且只能包含 A、C、G、T。
- **输出：** `consensus`、`max_score`、`min_score` 和 `scores` 列表，每项包含 `kmer` 与 `score`。
- **示例代码：**

```python
from dnakit.regulation import score_motif_kmers

result = score_motif_kmers(
    [[0, 0, 10, 0], [0, 0, 0, 10], [0, 0, 0, 10], [0, 0, 0, 10]],
    ["GTTT", "GTTA"],
)
print(result["consensus"], round(result["max_score"], 3))
for item in result["scores"]:
    print(item["kmer"], round(item["score"], 3))
```

- **示例输出：**

```text
GTTT 5.554
GTTT 5.554
GTTA -0.4
```

## 4) `dnakit.regulation.read_fragment_qc()`

- **作用：** 计算每个细胞条码的 fragment 数、重复比例和线粒体比例。
- **处理方法：** 读取纯文本、无表头的五列 paired-end fragment 文件（染色体、起点、终点、条码、读数），逐行校验后调用 snapATAC2 `pp.import_fragments()`，再读取其 `obs` 中的 `n_fragment`、`frac_dup` 和 `frac_mito`。`frac_dup` 是重复读数占全部读数的比例，`frac_mito` 是 chrM/M 上的 fragment 占全部 fragment 的比例。
- **API：** `dnakit.regulation.read_fragment_qc(fragment_file, chrom_sizes, *, min_num_fragments=1, max_fragments=5000000, max_cells=100000)`。
- **输入：** fragment 文件；`chrom_sizes` 为 `{染色体: 长度}`，文件中的染色体必须都在其中。`min_num_fragments` 会传给 snapATAC2，并去掉 fragment 数少于该值的条码。
- **输出：** `cell_count` 与 `cells` 列表，每项包含 `barcode`、`n_fragment`、`frac_dup` 和 `frac_mito`。
- **示例代码：**

```python
from dnakit.regulation import read_fragment_qc

result = read_fragment_qc(
    "fragments.tsv",
    {"chr1": 2000, "chr2": 2000, "chrM": 1000},
    min_num_fragments=1,
)
for cell in result["cells"]:
    print(cell["barcode"], cell["n_fragment"], round(cell["frac_dup"], 3), cell["frac_mito"])
```

- **示例输出：**

```text
cellA 3 0.2 0.25
cellB 3 0.4 0.0
```

## 5) `dnakit.regulation.count_tile_matrix()`

- **作用：** 生成细胞 × 固定宽度 bin 的计数矩阵。
- **处理方法：** 调用 snapATAC2 `pp.add_tile_matrix()`。bin 名称为 `chrom:start-end`（0-based，半开）。默认排除 chrM、chrY、M 和 Y。`counting_strategy` 取 `"fragment"`（每个 fragment 计一次）、`"insertion"` 或 `"paired-insertion"`。
- **API：** `dnakit.regulation.count_tile_matrix(fragment_file, chrom_sizes, *, bin_size=500, counting_strategy="fragment", exclude_chroms=("chrM", "chrY", "M", "Y"), min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_entries=1000000)`。
- **输入：** 与 `read_fragment_qc()` 相同的 fragment 文件和 `chrom_sizes`；`bin_size` 为正整数，单位为 bp。
- **输出：** `shape`、`barcodes`、`features`、`nonzero_count`、`total_count` 和 `entries`。`entries` 是按行、列排序的 `[细胞序号, bin 序号, 计数]` 三元组，序号从 0 开始。
- **示例代码：**

```python
from dnakit.regulation import count_tile_matrix

result = count_tile_matrix(
    "fragments.tsv",
    {"chr1": 2000, "chr2": 2000, "chrM": 1000},
    bin_size=500,
    min_num_fragments=1,
)
print(result["shape"], result["total_count"])
print(result["features"][:4])
print(result["entries"])
```

- **示例输出：**

```text
[2, 8] 6
['chr1:0-500', 'chr1:500-1000', 'chr1:1000-1500', 'chr1:1500-2000']
[[0, 0, 2], [0, 1, 1], [1, 1, 1], [1, 4, 1], [1, 6, 1]]
```

## 6) `dnakit.regulation.count_peak_fragments()`

- **作用：** 生成细胞 × 峰的计数矩阵。
- **处理方法：** 调用 snapATAC2 `pp.make_peak_matrix()`。峰来自 BED 文件的前三列（染色体、起点、终点），额外列被忽略。峰名称为 `chrom:start-end`。
- **API：** `dnakit.regulation.count_peak_fragments(fragment_file, chrom_sizes, peak_file, *, counting_strategy="fragment", min_num_fragments=1, max_fragments=5000000, max_cells=100000, max_peaks=100000, max_entries=1000000)`。
- **输入：** fragment 文件、`chrom_sizes` 和 BED 峰文件；BED 中的染色体也必须在 `chrom_sizes` 中。
- **输出：** 与 `count_tile_matrix()` 相同的字段，`features` 为峰名称。
- **示例代码：**

```python
from dnakit.regulation import count_peak_fragments

result = count_peak_fragments(
    "fragments.tsv",
    {"chr1": 2000, "chr2": 2000, "chrM": 1000},
    "peaks.bed",
    min_num_fragments=1,
)
print(result["features"])
print(result["entries"])
```

- **示例输出：**

```text
['chr1:50-400', 'chr1:590-950', 'chr2:1100-1500']
[[0, 0, 2], [0, 1, 1], [1, 1, 1], [1, 2, 1]]
```

## 7) `dnakit.regulation.predict_sequences()`

- **作用：** 用 tangermeme 的 `predict()` 对一组等长 DNA 序列做批量预测。
- **处理方法：** 把序列转换为 one-hot（A、C、G、T 的顺序），按 `batch_size` 调用模型。模型必须是 `torch.nn.Module` 且处于 eval 模式（否则报错，因为 dropout 与 batch norm 会改变输出）；只支持单个输出张量。
- **API：** `dnakit.regulation.predict_sequences(model, sequences, *, batch_size=32, max_sequences=256, max_sequence_length=10000)`。
- **输入：** PyTorch 模型和序列（格式与 `scan_motif_hits()` 相同）。所有序列必须长度相同，且只能包含 A、C、G、T。
- **输出：** `sequence_count`、`shape` 和 `predictions`（嵌套列表，第一维为序列）。
- **示例代码：**

```python
import torch
from dnakit.regulation import predict_sequences

torch.manual_seed(0)
model = torch.nn.Sequential(torch.nn.Flatten(), torch.nn.Linear(4 * 6, 1)).eval()
result = predict_sequences(model, {"a": "ACGTAC", "b": "TTTTGG"})
print(result["shape"])
print([round(row[0], 4) for row in result["predictions"]])
```

- **示例输出：**

```text
[2, 1]
[0.094, 0.1866]
```

## 8) `dnakit.regulation.marginalize_motif_effect()`

- **作用：** 边缘化分析：比较把同一个基序写入每条背景序列前后，模型输出的变化。
- **处理方法：** 调用 tangermeme 的 `marginalize()`。基序覆盖到每条序列上；`start` 为空时，基序放在序列中央，位置为 `长度 // 2 - 基序长度 // 2`。`delta` 是写入后的输出减去写入前的输出。
- **API：** `dnakit.regulation.marginalize_motif_effect(model, sequences, motif, *, start=None, max_sequences=256, max_sequence_length=10000)`。
- **输入：** 模型、序列和一个 ACGT 基序字符串。`start` 是基序在每条序列中的 0-based 起点，必须落在序列内。
- **输出：** `motif`、`start`、`shape`、`y_before`、`y_after` 和 `delta`。
- **示例代码：**

```python
import torch
from dnakit.regulation import marginalize_motif_effect


class CountA(torch.nn.Module):
    """Outputs the number of A bases in each sequence (one-hot channel 0)."""

    def forward(self, inputs):
        return inputs[:, 0, :].sum(dim=1, keepdim=True)


model = CountA().eval()
result = marginalize_motif_effect(model, {"c": "CCCCCCCCCC", "ac": "ACCCCCCCCC"}, "AAAA")
print(result["start"], result["y_before"], result["y_after"], result["delta"])
```

- **示例输出：**

```text
3 [[0.0], [1.0]] [[4.0], [5.0]] [[4.0], [4.0]]
```

## 9) `dnakit.regulation.deep_lift_shap_attribution()`

- **作用：** 计算模型某个输出对每个位置、每个碱基的 DeepLIFT/SHAP 归因。
- **处理方法：** 调用 tangermeme 的 `deep_lift_shap()`。参考序列是每条输入序列的双核苷酸打乱版本，由 `random_state` 固定，因此结果可以复现。`max_positions` 限制序列数乘以长度。
- **API：** `dnakit.regulation.deep_lift_shap_attribution(model, sequences, *, target=0, n_shuffles=20, random_state=0, max_sequences=256, max_sequence_length=10000, max_positions=200000)`。
- **输入：** 模型（eval 模式，单个输出张量）、序列，以及要解释的输出列 `target`。
- **输出：** `shape` 和 `attributions`，按序列、位置、碱基（A、C、G、T）排列。
- **示例代码：**

```python
import torch
from dnakit.regulation import deep_lift_shap_attribution


class CountA(torch.nn.Module):
    """Outputs the number of A bases in each sequence (one-hot channel 0)."""

    def forward(self, inputs):
        return inputs[:, 0, :].sum(dim=1, keepdim=True)


model = CountA().eval()
result = deep_lift_shap_attribution(model, {"s": "ACGTTGCAAC"}, n_shuffles=4, random_state=3)
print(result["shape"])
print([round(value, 2) for value in result["attributions"][0][0]])
```

- **示例输出：**

```text
[1, 10, 4]
[0.0, -0.0, -0.0, -0.0]
```

## 安装可选依赖

```bash
python -m pip install "dnakit[bio-gimmemotifs]"     # gimmemotifs，read_motif_file / scan_motif_hits / score_motif_kmers
python -m pip install "dnakit[bio-snapatac2]"       # snapATAC2，read_fragment_qc / count_tile_matrix / count_peak_fragments
python -m pip install "dnakit[bio-tangermeme]"      # tangermeme 与 PyTorch，predict_sequences / marginalize_motif_effect / deep_lift_shap_attribution
```

各包的项目地址：[gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs)、[snapATAC2](https://github.com/scverse/snapatac2)、[tangermeme](https://github.com/jmschrei/tangermeme)。

## 验证方式

测试文件是 `others/tests/unit/test_regulation_backends.py`，共 20 项，主要检查：

- gimmemotifs：PFM、TRANSFAC、JASPAR 三种格式读取出相同的计数；反向互补链的命中坐标与手工计算一致；对数优势得分与 `ln(ppm / 0.25 + 0.01)` 的手工求和以及直接调用 `Motif.score_kmer()` 一致；超过 `max_hits_per_sequence` 时报错。
- snapatac2：fragment QC（细胞 fragment 数、重复比例、线粒体比例）与手算值一致，并与直接调用 `pp.import_fragments()` 的结果一致；bin 计数矩阵与手工计数一致，并与直接调用 `add_tile_matrix()` 的结果一致；峰计数矩阵与手工计数一致。
- tangermeme：预测结果与直接调用模型一致；基序边缘化的前后输出差值与手工计数一致；DeepLIFT 归因与直接调用 `deep_lift_shap()` 一致，并且相同种子下可复现。
- 文件校验、分隔符、染色体名称、坐标和行数上限都会被拒绝；缺少依赖时给出安装提示；上游抛出的 `BaseException`（例如 pyo3 panic）会转换为 `BackendExecutionError`（测试中用替身异常验证）。

## 当前范围

这一批包装三个包中的 DNA 层面操作。已包装的是：gimmemotifs 的 PFM/TRANSFAC/JASPAR 读取、基于 `scan_all()` 的序列扫描和 k-mer 打分；snapATAC2 的 fragment QC、bin 计数矩阵和峰计数矩阵；tangermeme 的批量预测、基序边缘化和 DeepLIFT/SHAP 归因。

未包装的操作及原因如下：

- gimmemotifs 的 `Scanner` 与基于 FPR 的阈值：需要基因组或背景序列文件，本批只做离线扫描。
- gimmemotifs 的基序发现程序（MEME、HOMER 等）：依赖外部可执行程序；基序写出函数也未包装。
- snapatac2 的单端 fragment 格式（`is_paired=False`）：需要额外的链信息列，本批只接受五列 paired-end 格式；gzip 输入、基因活性矩阵、TSS 富集（需要基因注释）、peak calling、降维和聚类也未包装。
- tangermeme 的饱和突变、变异效应、ablation、序列设计以及 `read_meme` 等读取函数：可在后续批次中按需添加。
