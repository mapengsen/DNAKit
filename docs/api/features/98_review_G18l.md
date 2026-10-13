# G18l：cooler、cooltools 与 pairtools 的评审函数

这一组对 cooler、cooltools 和 pairtools 中原先状态为 `other` 或 `cap` 的函数逐一评审，共包装 24 个 DNA 或基因组层面、可离线运行的函数：cooler 4 个、cooltools 16 个、pairtools 4 个。其余函数按原因记录在账本中，例如 `cli`、`internal`、`not-dna`、`needs-network`、`family-variant` 或 `duplicate-of`。
函数位于 `dnakit.epigenome` 的三个评审模块中（`cooler_review_functions.py`、`cooltools_review_functions.py`、`pairtools_review_functions.py`）。每个函数只在调用时导入上游包；缺少可选依赖时抛出 `BackendUnavailableError`。输入先检查，上游的 `ValueError`、`KeyError`、`OSError` 等转为 `ConfigurationError`。所有公开函数都有关键字参数 `verbose=True`，返回普通字典并附 `provenance`。
特征分解与诱饵函数（`dense_trans_eigenvectors`、`replace_cis_with_decoys`）在调用期间保存并恢复 numpy 全局随机状态与 C 库 `random()` 的状态（glibc 与 BSD 提供 `initstate`/`setstate`），并为 ARPACK 固定初始向量；这一补丁在 `finally` 中恢复，但不是线程安全的。在没有 `initstate`/`setstate` 的平台上，C 随机状态不会被隔离。`call_chromatin_dots` 在 100 bp 等细分辨率下必须给出核宽，因此示例的核宽显式设定。
示例都使用同一个合成 cooler（`synthetic.cool`，两条染色体，长 2000 bp 与 1500 bp，分辨率 100 bp，第一个示例准备步骤生成），并在当前目录下写出输出文件。示例输出由运行代码得到。
## 示例文件

以下文件在运行示例前写入示例目录。

**示例文件：** `track.bedgraph`

```text
chr1	0	100	1.0
chr1	100	200	2.0
chr1	200	300	3.0
chr2	0	100	0.5
chr2	100	200	0.6
```

**示例文件：** `features.bed`

```text
chr1	300	500
chr2	100	200
```

**示例文件：** `features.bedpe`

```text
chr1	100	150	chr1	700	750
chr1	500	550	chr2	50	80
```

**示例文件：** `view.bed`

```text
chr1	0	400	A
chr1	600	1000	B
chr2	0	800	X
```

**示例文件：** `cells.pairs`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	100	chr1	300	+	-	UU
r2	chr1	100	chr1	300	+	-	UU
r3	chr2	10	chr2	20	+	-	UU
r4	chr1	500	chr2	100	+	-	UU
```

**示例文件：** `cov.pairs`

```text
## pairs format v1.0
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	100	chr1	300	+	-	UU
r2	chr1	100	chr1	300	+	-	UU
r3	chr1	100	chr1	300	+	-	UU
r4	chr2	10	chr2	20	+	-	UU
r5	chr1	500	chr2	100	+	-	UU
```

**示例文件：** `reads.sam`

```text
@HD	VN:1.6	SO:queryname
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:800
r1	65	chr1	101	60	10M	=	101	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r1	129	chr2	301	60	10M	=	301	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r2	65	chr2	51	60	10M	=	51	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
r2	129	chr2	201	60	10M	=	201	0	ACGTACGTAC	IIIIIIIIII	NM:i:0	AS:i:10	XS:i:0
```


## 示例准备

这个示例创建 `synthetic.cool`，并删除以前运行留下的输出文件。

**示例代码：**

```python
import os

import cooler
import numpy as np
import pandas as pd

for name in ["synthetic.cool", "cells.scool", "all.pairs", "strict.pairs", "removed.pairs",
             "reads.pairs"]:
    if os.path.exists(name):
        os.remove(name)
rng = np.random.default_rng(0)
bins = cooler.binnify(pd.Series({"chr1": 2000, "chr2": 1500}), 100)
rows = []
for i in range(len(bins)):
    for j in range(i, len(bins)):
        same = bins.chrom[i] == bins.chrom[j]
        count = rng.poisson(6.0 / (1 + abs(i - j)) if same else 0.5)
        if i == j:
            count += 3
        if count > 0:
            rows.append((i, j, int(count)))
rows.append((2, 9, 40))
pixels = pd.DataFrame(rows, columns=["bin1_id", "bin2_id", "count"])
pixels = pixels.sort_values(["bin1_id", "bin2_id"])
cooler.create_cooler("synthetic.cool", bins, pixels, dtypes={"count": "int"}, ensure_sorted=True)
cooler.balance_cooler(cooler.Cooler("synthetic.cool"), store=True)
clr = cooler.Cooler("synthetic.cool")
print("bins", len(clr.bins()[:]), "chromosomes", list(clr.chromnames))
```

**示例输出：**

```text
bins 35 chromosomes ['chr1', 'chr2']
```

## cooler 函数（4 个）

## 1) `dnakit.epigenome.create_single_cell_cooler()`

- **作用：** 把每个细胞的 .pairs 文件装入一个单细胞 cooler（.scool）文件。
- **处理方法：** 每个细胞的 .pairs 按 binsize 分箱，上三角像素按 bin 对求和，然后调用 cooler.create_scool。输出为新文件，不覆盖已有文件。
- **API：** `dnakit.epigenome.create_single_cell_cooler(output_path, pairs_by_cell, chromsizes, binsize, *, verbose=True)`.
- **输入：** `output_path`：新的 .scool 文件；`pairs_by_cell`：细胞名到 .pairs 文件的映射；`chromsizes`：染色体长度；`binsize`：分箱大小。
- **输出：** 含 `bin_count`、`cell_count` 和 `cells`（每个细胞的 `pair_count`、`nonzero_pixels`）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import create_single_cell_cooler

result = create_single_cell_cooler(
    "cells.scool",
    {"cellA": "cells.pairs", "cellB": "cells.pairs"},
    {"chr1": 1000, "chr2": 800},
    100,
    verbose=False,
)
print(result["operation"], result["bin_count"], result["cell_count"])
for cell in result["cells"]:
    print(cell["name"], cell["pair_count"], cell["nonzero_pixels"])
```

- **示例输出：**

```text
create_single_cell_cooler 18 2
cellA 4 3
cellB 4 3
```

## 2) `dnakit.epigenome.list_cooler_groups()`

- **作用：** 列出 HDF5 文件中保存 cooler 数据集的组路径。
- **处理方法：** 调用 cooler.api.list_coolers，结果按组路径排序。普通 .cool 文件返回 `/`，多分辨率文件返回各分辨率组。
- **API：** `dnakit.epigenome.list_cooler_groups(path, *, verbose=True)`.
- **输入：** `path`：HDF5 文件（.cool、.mcool 或 .scool）。
- **输出：** 含 `group_count` 和 `groups` 的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import list_cooler_groups

print(list_cooler_groups("synthetic.cool", verbose=False)["groups"])
print(list_cooler_groups("cells.scool", verbose=False)["groups"])
```

- **示例输出：**

```text
['/']
['/cells/cellA', '/cells/cellB']
```

## 3) `dnakit.epigenome.list_scool_cells()`

- **作用：** 列出单细胞（.scool）文件中的细胞。
- **处理方法：** 先用 cooler.fileops.is_scool_file 检查文件类型，非单细胞文件报错；然后调用 list_scool_cells，按组路径排序，并给出细胞名。
- **API：** `dnakit.epigenome.list_scool_cells(path, *, verbose=True)`.
- **输入：** `path`：单细胞 cooler 文件。
- **输出：** 含 `cell_count` 和 `cells`（每项为 `name` 与 `group`）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import list_scool_cells

result = list_scool_cells("cells.scool", verbose=False)
print(result["cell_count"], [cell["name"] for cell in result["cells"]])
```

- **示例输出：**

```text
2 ['cellA', 'cellB']
```

## 4) `dnakit.epigenome.infer_cooler_binsize()`

- **作用：** 从 bin 表推断 cooler 文件的分箱大小。
- **处理方法：** 调用 cooler.util.get_binsize：除每条染色体的最后一个 bin 外，其余 bin 必须等长；不等长时 `binsize` 为 None。
- **API：** `dnakit.epigenome.infer_cooler_binsize(cool_path, *, verbose=True)`.
- **输入：** `cool_path`：cooler 文件。
- **输出：** 含 `binsize`、`uniform` 和 `bin_count` 的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import infer_cooler_binsize

result = infer_cooler_binsize("synthetic.cool", verbose=False)
print(result["binsize"], result["uniform"], result["bin_count"])
```

- **示例输出：**

```text
100 True 35
```

## cooltools 函数（16 个）

## 5) `dnakit.epigenome.call_chromatin_dots()`

- **作用：** 用指定宽度的 HiCCUPS 型核（donut、vertical、horizontal、lowleft）调用染色质环（dots）。
- **处理方法：** 用 cooltools.expected_cis 计算顺式期望（balanced.avg），用 numutils.get_kernel 生成四个核，再调用 cooltools.dots。100 bp 等细分辨率不再推荐默认核，因此需要显式给出核宽。
- **API：** `dnakit.epigenome.call_chromatin_dots(cool_path, *, weight='weight', kernel_half_width=3, donut_inner_half_width=1, max_loci_separation=10000000, n_lambda_bins=40, verbose=True)`.
- **输入：** `cool_path`：已平衡的 cooler；`kernel_half_width`、`donut_inner_half_width`：核的外半宽与内半宽。
- **输出：** 含 `dot_count` 和 `dots`（坐标、计数与各核的 la_exp 值与 qval）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import call_chromatin_dots

result = call_chromatin_dots("synthetic.cool", verbose=False)
print(result["kernel_half_width"], result["donut_inner_half_width"], result["dot_count"])
```

- **示例输出：**

```text
3 1 0
```

## 6) `dnakit.epigenome.align_track_to_cooler_bins()`

- **作用：** 把 bedGraph 式的 track 对齐到 cooler 的 bin 上。
- **处理方法：** 调用 cooltools.lib.common.align_track_with_cooler：检查 bin 大小，合并 track 与 bin 表；权重为非有限值的坏 bin 设为 None（可关闭）。
- **API：** `dnakit.epigenome.align_track_to_cooler_bins(track_path, cool_path, *, weight='weight', mask_bad_bins=True, verbose=True)`.
- **输入：** `track_path`：含 chrom、start、end、value 四列的 TSV；`cool_path`：已平衡的 cooler。
- **输出：** 含 `bin_count`、`valued_bin_count` 和 `bins` 的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import align_track_to_cooler_bins

result = align_track_to_cooler_bins("track.bedgraph", "synthetic.cool", verbose=False)
print(result["bin_count"], result["valued_bin_count"])
for row in result["bins"]:
    if row["value"] is not None:
        print(row["chrom"], row["start"], row["end"], row["value"])
```

- **示例输出：**

```text
35 5
chr1 0 100 1.0
chr1 100 200 2.0
chr1 200 300 3.0
chr2 0 100 0.5
chr2 100 200 0.6
```

## 7) `dnakit.epigenome.make_track_view()`

- **作用：** 从 track 生成视窗（viewframe），每条染色体一个区域。
- **处理方法：** 调用 cooltools.lib.common.view_from_track：每条染色体取最小起点与最大终点，并用 bioframe 生成区域名。
- **API：** `dnakit.epigenome.make_track_view(track_path, *, verbose=True)`.
- **输入：** `track_path`：含 chrom、start、end 的 TSV。
- **输出：** 含 `region_count` 和 `regions`（chrom、start、end、name）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import make_track_view

result = make_track_view("track.bedgraph", verbose=False)
for row in result["regions"]:
    print(row["chrom"], row["start"], row["end"], row["name"])
```

- **示例输出：**

```text
chr1 0 300 chr1
chr2 0 200 chr2
```

## 8) `dnakit.epigenome.mask_track_bad_bins()`

- **作用：** 把落在 cooler 坏 bin 上的 track 值设为缺失。
- **处理方法：** 调用 cooltools.lib.common.mask_cooler_bad_bins：只保留 chrom、start、end 与 bin 表完全一致的行，权重非有限的 bin 其值设为 NaN。
- **API：** `dnakit.epigenome.mask_track_bad_bins(track_path, cool_path, *, weight='weight', verbose=True)`.
- **输入：** `track_path`：TSV；`cool_path`：cooler 文件；`weight`：权重列名。
- **输出：** 含 `row_count`、`masked_count` 和 `rows` 的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import mask_track_bad_bins

result = mask_track_bad_bins("track.bedgraph", "synthetic.cool", verbose=False)
print(result["row_count"], result["masked_count"])
```

- **示例输出：**

```text
5 0
```

## 9) `dnakit.epigenome.assign_paired_features_to_view()`

- **作用：** 为成对特征（BEDPE）的两端分别标注视窗区域名。
- **处理方法：** 读取六列 BEDPE 与 BED 视窗，调用 cooltools.lib.common.assign_view_paired；落在任何区域之外的端点为 None。
- **API：** `dnakit.epigenome.assign_paired_features_to_view(features_path, view_path, *, drop_unassigned=False, verbose=True)`.
- **输入：** `features_path`：六列 BEDPE；`view_path`：视窗 BED（可有名字列）；`drop_unassigned`：是否丢弃未标注的特征。
- **输出：** 含 `feature_count`、`assigned_count` 和 `features`（含 `region1`、`region2`）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import assign_paired_features_to_view

result = assign_paired_features_to_view("features.bedpe", "view.bed", verbose=False)
for row in result["features"]:
    print(row["chrom1"], row["start1"], row["chrom2"], row["start2"],
          row["region1"], row["region2"])
```

- **示例输出：**

```text
chr1 100 chr1 700 A B
chr1 500 chr2 50 None X
```

## 10) `dnakit.epigenome.assign_bins_to_regions()`

- **作用：** 给出每个 bin 所在的区域名（按 bin 下标的跨度判断）。
- **处理方法：** 按 bin_start 排序区域，拒绝重叠区域，然后调用 cooltools.lib.common.assign_regions_to_bins；不在任何区域内的 bin 为 None。
- **API：** `dnakit.epigenome.assign_bins_to_regions(bin_ids, regions, *, verbose=True)`.
- **输入：** `bin_ids`：bin 下标列表；`regions`：每项含 name、bin_start、bin_end 的映射列表。
- **输出：** 含 `bin_count`、`assigned_count` 和 `assignments`（bin_id 与 region）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import assign_bins_to_regions

result = assign_bins_to_regions(
    [0, 4, 5, 11, 12],
    [
        {"name": "a", "bin_start": 0, "bin_end": 5},
        {"name": "b", "bin_start": 5, "bin_end": 12},
    ],
    verbose=False,
)
print([(row["bin_id"], row["region"]) for row in result["assignments"]])
```

- **示例输出：**

```text
[(0, 'a'), (4, 'a'), (5, 'b'), (11, 'b'), (12, None)]
```

## 11) `dnakit.epigenome.dense_trans_eigenvectors()`

- **作用：** 计算稠密全基因组矩阵的反式（trans）区室特征向量。
- **处理方法：** 调用 cooltools.api.eigdecomp.trans_eig：去除顺式块，截断极端值，用顺式伪数据（fake_cis）补全，再做特征分解；特征向量每列最大绝对值为正。调用期间 numpy 全局随机状态与 C 库 random() 的状态被保存并恢复，ARPACK 使用固定初始向量。
- **API：** `dnakit.epigenome.dense_trans_eigenvectors(matrix, partition, *, n_eigs=3, perc_top=99.95, perc_bottom=1.0, verbose=True)`.
- **输入：** `matrix`：对称的方阵；`partition`：从 0 递增到矩阵大小的染色体偏移；`n_eigs`：特征向量个数。
- **输出：** 含 `eigenvalues` 和 `eigenvectors` 的字典，附 `provenance`。
- **示例代码：**

```python
import numpy as np

from dnakit.epigenome import dense_trans_eigenvectors

rng = np.random.default_rng(5)
base = rng.random((12, 12))
matrix = (base + base.T).tolist()
result = dense_trans_eigenvectors(matrix, [0, 6, 12], n_eigs=2, verbose=False)
print(result["size"], result["partition"])
print([round(value, 4) for value in result["eigenvalues"]])
```

- **示例输出：**

```text
12 [0, 6, 12]
[2.5948, -2.4452]
```

## 12) `dnakit.epigenome.pairwise_diagonal_sums()`

- **作用：** 计算每对染色体臂内部各对角线的像素和（顺式块）。
- **处理方法：** 调用 cooltools.api.expected.diagsum_pairwise，按区域对与距离排序；默认忽略最近的 2 条对角线。
- **API：** `dnakit.epigenome.pairwise_diagonal_sums(cool_path, *, weight='weight', ignore_diags=2, verbose=True)`.
- **输入：** `cool_path`：已平衡的 cooler；`weight`、`ignore_diags`。
- **输出：** 含 `row_count` 和 `rows`（region1、region2、dist、n_total、n_valid、count.sum）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import pairwise_diagonal_sums

result = pairwise_diagonal_sums("synthetic.cool", verbose=False)
chr1_rows = [row for row in result["rows"] if row["region1"] == "chr1"]
print(result["row_count"], len(chr1_rows))
```

- **示例输出：**

```text
35 20
```

## 13) `dnakit.epigenome.pairwise_block_sums()`

- **作用：** 计算每对染色体之间矩形块的像素和（反式块）。
- **处理方法：** 调用 cooltools.api.expected.blocksum_pairwise，结果按区域对排序。
- **API：** `dnakit.epigenome.pairwise_block_sums(cool_path, *, weight='weight', verbose=True)`.
- **输入：** `cool_path`：已平衡的 cooler；`weight`。
- **输出：** 含 `row_count` 和 `rows`（region1、region2、n_valid、count.sum）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import pairwise_block_sums

result = pairwise_block_sums("synthetic.cool", verbose=False)
row = result["rows"][0]
print(result["row_count"], row["region1"], row["region2"], row["n_valid"], row["count.sum"])
```

- **示例输出：**

```text
1 chr1 chr2 280 134
```

## 14) `dnakit.epigenome.smooth_contact_decay()`

- **作用：** 在 log 空间中平滑平衡后的接触-距离曲线。
- **处理方法：** 用 cooltools.expected_cis 计算顺式接触-距离表，再调用 genomewide_smooth_cvd（mode="genomewide"，全部区域合并）或 per_region_smooth_cvd（mode="per_region"，逐区域）。后缀固定为 .smoothed。
- **API：** `dnakit.epigenome.smooth_contact_decay(cool_path, *, mode='genomewide', weight='weight', sigma_log10=0.1, window_sigma=5, points_per_sigma=10, verbose=True)`.
- **输入：** `cool_path`：已平衡的 cooler；`mode`；`sigma_log10`、`window_sigma`、`points_per_sigma`：高斯核参数。
- **输出：** 含 `row_count` 和 `rows`（区域、距离与平滑列）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import smooth_contact_decay

result = smooth_contact_decay("synthetic.cool", verbose=False)
row = result["rows"][2]
print(result["mode"], result["row_count"], row["dist"], round(row["balanced.avg.smoothed"], 6))
```

- **示例输出：**

```text
genomewide 35 2 0.07744
```

## 15) `dnakit.epigenome.saddle_compartment_strength()`

- **作用：** 从 saddle 求和矩阵与计数矩阵计算区室强度。
- **处理方法：** 调用 cooltools.api.saddle.saddle_strength：对 k 个 bin 的两端角，计算同区室（AA+BB）与异区室（AB+BA）均值之比。
- **API：** `dnakit.epigenome.saddle_compartment_strength(saddle_sums, saddle_counts, *, verbose=True)`.
- **输入：** `saddle_sums`、`saddle_counts`：同形状的方阵（可含 None）。
- **输出：** 含 `size` 和 `strengths`（每项为 k 与 strength）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import saddle_compartment_strength

sums = [[1, 2, 3, 4], [2, 5, 6, 7], [3, 6, 8, 9], [4, 7, 9, 10]]
counts = [[2, 2, 2, 2] for _ in range(4)]
result = saddle_compartment_strength(sums, counts, verbose=False)
print([(row["k"], row["strength"]) for row in result["strengths"]])
```

- **示例输出：**

```text
[(1, 1.375), (2, 1.15), (3, 1.03)]
```

## 16) `dnakit.epigenome.expand_features_to_bins()`

- **作用：** 把 BED 或 BEDPE 特征按侧翼长度扩展，并对齐到 bin 边界。
- **处理方法：** 调用 cooltools.api.snipping.expand_align_features：中心点为 floor(中心/分辨率)，两侧各加 flank//分辨率 个 bin，边界乘以分辨率。
- **API：** `dnakit.epigenome.expand_features_to_bins(features_path, *, flank, resolution, format='bed', verbose=True)`.
- **输入：** `features_path`：BED 或 BEDPE；`flank`：侧翼长度（bp）；`resolution`：分辨率（bp）；`format`：`bed` 或 `bedpe`。
- **输出：** 含 `feature_count` 和 `features`（含原始与扩展后的坐标）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import expand_features_to_bins

result = expand_features_to_bins("features.bed", flank=100, resolution=50, verbose=False)
for row in result["features"]:
    print(row["chrom"], row["start"], row["end"])
```

- **示例输出：**

```text
chr1 300 550
chr2 50 300
```

## 17) `dnakit.epigenome.adaptive_coarsegrain_matrix()`

- **作用：** 自适应地粗粒化接触矩阵：稀疏像素合并到邻域中直到计数足够。
- **处理方法：** 调用 cooltools.lib.numutils.adaptive_coarsegrain：每次做 2×2 合并，直到每个像素的原始计数不少于 cutoff，或达到 max_levels 层。输出与输入同形状。
- **API：** `dnakit.epigenome.adaptive_coarsegrain_matrix(matrix, counts, *, cutoff=5.0, max_levels=8, min_shape=8, verbose=True)`.
- **输入：** `matrix`：接触矩阵；`counts`：同形状的原始计数矩阵；`cutoff`、`max_levels`、`min_shape`。
- **输出：** 含 `size` 和 `matrix` 的字典，附 `provenance`。
- **示例代码：**

```python
import numpy as np

from dnakit.epigenome import adaptive_coarsegrain_matrix

rng = np.random.default_rng(3)
matrix = rng.random((9, 9))
counts = np.full((9, 9), 2.0)
result = adaptive_coarsegrain_matrix(
    matrix.tolist(), counts.tolist(), cutoff=3, max_levels=2, min_shape=4, verbose=False
)
print(result["size"], len(result["matrix"]), len(result["matrix"][0]))
```

- **示例输出：**

```text
9 9 9
```

## 18) `dnakit.epigenome.replace_cis_with_decoys()`

- **作用：** 把掩码标记的顺式像素替换为伪造的诱饵接触（decoy）。
- **处理方法：** 调用 cooltools.lib.numutils.fake_cis：掩码为 1 的像素用 C 库 random() 生成的值替换，掩码为 0 的像素保持不变。调用前将 C 随机流设为由 seed 决定的私有状态，调用后恢复。
- **API：** `dnakit.epigenome.replace_cis_with_decoys(matrix, mask, *, seed=0, verbose=True)`.
- **输入：** `matrix`：方阵；`mask`：同形状的 0/1 矩阵；`seed`：随机种子。
- **输出：** 含 `replaced_pixels`、`seed` 和 `matrix` 的字典，附 `provenance`。
- **示例代码：**

```python
import numpy as np

from dnakit.epigenome import replace_cis_with_decoys

rng = np.random.default_rng(4)
base = rng.random((12, 12))
matrix = (base + base.T).tolist()
mask = [[1 if i < 6 and j < 6 else 0 for j in range(12)] for i in range(12)]
result = replace_cis_with_decoys(matrix, mask, seed=0, verbose=False)
print(result["replaced_pixels"], result["seed"])
print(result["matrix"][6][6] == matrix[6][6], result["matrix"][0][0] == matrix[0][0])
```

- **示例输出：**

```text
36 0
True False
```

## 19) `dnakit.epigenome.expected_full_table()`

- **作用：** 计算每一对染色体区域的期望接触数（顺式为距离衰减，反式为块均值）。
- **处理方法：** 调用 cooltools 沙盒中的 expected_full（fast=True 时为 expected_full_fast）。输出按区域对与距离排序。
- **API：** `dnakit.epigenome.expected_full_table(cool_path, *, weight='weight', fast=False, verbose=True)`.
- **输入：** `cool_path`：已平衡的 cooler；`weight`；`fast`：是否使用较快的变体。
- **输出：** 含 `row_count` 和 `rows`（region1、region2、dist、n_valid、count.sum、expected）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import expected_full_table

result = expected_full_table("synthetic.cool", verbose=False)
print(result["row_count"], result["rows"][0]["region1"], result["rows"][0]["region2"])
```

- **示例输出：**

```text
36 chr1 chr1
```

## 20) `dnakit.epigenome.observed_over_expected_pixels()`

- **作用：** 计算 cooler 每个像素的观测/期望值。
- **处理方法：** 用 expected_full 计算期望，将 bin 标注到视窗区域，再调用沙盒的 obs_over_exp：oe = 平衡计数 / 期望值。
- **API：** `dnakit.epigenome.observed_over_expected_pixels(cool_path, *, weight='weight', verbose=True)`.
- **输入：** `cool_path`：已平衡的 cooler；`weight`。
- **输出：** 含 `pixel_count` 和 `pixels`（bin 下标、平衡计数、期望与 oe）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import observed_over_expected_pixels

result = observed_over_expected_pixels("synthetic.cool", verbose=False)
valid = [row for row in result["pixels"] if row["oe"] is not None]
print(result["pixel_count"], len(valid))
print(round(valid[0]["oe"], 4))
```

- **示例输出：**

```text
337 272
0.3723
```

## pairtools 函数（4 个）

## 21) `dnakit.epigenome.sam_to_pairs()`

- **作用：** 把每条 read 的 SAM 比对分类为 Hi-C 对，并写出 .pairs 文件。
- **处理方法：** 调用 pairtools.lib.parse.streaming_classify（常规解析器）；只保留唯一比对，写入标准 .pairs 表头（不含 SAM 文本）。输出为新文件。
- **API：** `dnakit.epigenome.sam_to_pairs(sam_path, output_path, *, min_mapq=40, max_molecule_size=500, max_inter_align_gap=20, walks_policy='mask', verbose=True)`.
- **输入：** `sam_path`：带 @SQ 头的 SAM 或 BAM；`output_path`：新的 .pairs 文件；`min_mapq`、`max_molecule_size`、`max_inter_align_gap`、`walks_policy`。
- **输出：** 含 `pair_count`、`chromosome_count` 和 `columns` 的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import sam_to_pairs

result = sam_to_pairs("reads.sam", "reads.pairs", verbose=False)
print(result["pair_count"], result["chromosome_count"])
print(open("reads.pairs").read())
```

- **示例输出：**

```text
2 2
## pairs format v1.0.0
#shape: upper triangle
#genome_assembly: unknown
#chromsize: chr1 1000
#chromsize: chr2 800
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
r1	chr1	101	chr2	301	+	+	UU
r2	chr2	51	chr2	201	+	+	UU
```

## 22) `dnakit.epigenome.filter_pairs_by_coverage()`

- **作用：** 删除片段端点被附近过多 pair 覆盖的 pair（高频相互作用端点）。
- **处理方法：** 调用 pairtools.lib.filterbycov.streaming_filterbycov：每个端点统计 max_dist 内的伙伴数，按 method 为 sum 或 max 汇总；超过 max_cov 的 pair 被删除，可写入 removed_path。
- **API：** `dnakit.epigenome.filter_pairs_by_coverage(pairs_path, output_path, *, method='max', max_dist=500, max_cov=8, removed_path=None, verbose=True)`.
- **输入：** `pairs_path`：.pairs 文件；`output_path`：保留的 pair 输出文件；`method`、`max_dist`（默认 500）、`max_cov`（默认 8）。
- **输出：** 含 `kept_count`、`removed_count` 的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import filter_pairs_by_coverage

everything = filter_pairs_by_coverage("cov.pairs", "all.pairs", max_cov=100, verbose=False)
print(everything["kept_count"], everything["removed_count"])
strict = filter_pairs_by_coverage(
    "cov.pairs",
    "strict.pairs",
    method="max",
    max_dist=500,
    max_cov=3,
    removed_path="removed.pairs",
    verbose=False,
)
print(strict["kept_count"], strict["removed_count"])
```

- **示例输出：**

```text
5 None
1 4
```

## 23) `dnakit.epigenome.mark_sam_duplicate()`

- **作用：** 给一条 SAM 记录设置重复标记（FLAG 位 1024），并把 Yt:Z 标签改为 DD。
- **处理方法：** 调用 pairtools.lib.dedup.mark_sam_as_dup（先把制表符换成 pairsam 字段分隔符 \x19，处理后再换回）。
- **API：** `dnakit.epigenome.mark_sam_duplicate(sam_line, *, verbose=True)`.
- **输入：** `sam_line`：一条至少 11 列的 SAM 文本记录。
- **输出：** 含 `sam`、`flag_before`、`flag_after` 与 `duplicate_flag_set` 的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import mark_sam_duplicate

line = "r1\t99\tchr1\t101\t60\t10M\t=\t201\t110\tACGTACGTAC\tIIIIIIIIII\tYt:Z:UU"
result = mark_sam_duplicate(line, verbose=False)
print(result["flag_before"], result["flag_after"], result["duplicate_flag_set"])
print(result["sam"])
```

- **示例输出：**

```text
99 1123 True
r1	1123	chr1	101	60	10M	=	201	110	ACGTACGTAC	IIIIIIIIII	Yt:Z:DD
```

## 24) `dnakit.epigenome.chromsizes_from_sam_header()`

- **作用：** 读取 SAM 或 BAM 头部 @SQ 行中的序列名与长度。
- **处理方法：** 用 pairtools 的 pysam 读取器打开文件，调用 headerops.get_chromsizes_from_pysam_header；保持头部顺序。
- **API：** `dnakit.epigenome.chromsizes_from_sam_header(sam_path, *, verbose=True)`.
- **输入：** `sam_path`：带 @SQ 头的 SAM 或 BAM 文件。
- **输出：** 含 `chromosome_count` 和 `chromosomes`（name、length）的字典，附 `provenance`。
- **示例代码：**

```python
from dnakit.epigenome import chromsizes_from_sam_header

result = chromsizes_from_sam_header("reads.sam", verbose=False)
for row in result["chromosomes"]:
    print(row["name"], row["length"])
```

- **示例输出：**

```text
chr1 1000
chr2 800
```
