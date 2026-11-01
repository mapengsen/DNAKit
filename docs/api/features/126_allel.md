# scikit-allel：群体遗传绘图与基因型存储

这一组函数是 `dnakit.population` 中对 [scikit-allel](https://github.com/cggh/scikit-allel) 的封装，覆盖群体遗传图形、纯合片段（ROH）和 VCF 转 Zarr。其中 15 个函数直接调用 scikit-allel 自带的绘图函数，绘制位点频谱、单倍型涂色与频率、两两距离与连锁不平衡矩阵以及变异定位图。它们写出一个新的 PNG 文件（150 dpi，Agg 后端，已存在的文件不会被覆盖），并返回路径和实际绘制的数值，不返回绘图对象。`roh_poissonhmm()` 调用 scikit-allel 的泊松隐马尔可夫模型来调用纯合片段，`allel_vcf_to_zarr()` 借助 [zarr](https://github.com/zarr-developers/zarr-python) 把 VCF 文件的数组写入 Zarr 存储。绘图函数只绘制传入的数据：它们自己不计算频谱、LD 或距离。

scikit-allel 使用 MIT 许可（已核对已安装包的元数据）。依赖都是可选的，只有调用对应函数时才会导入；缺少包时会报错并给出安装命令。示例从 `dnakit.population` 导入。

## 示例数据

只有最后一个函数读取输入文件：一个含 3 个变异、3 个样本的 VCF。其余示例直接传入小列表，并把 PNG 文件写入临时目录。

**示例文件：** `mini.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	v1	A	G	50	PASS	.	GT	0/0	0/1	1/1
chr1	250	v2	C	T	60	PASS	.	GT	0/1	0/0	0/1
chr1	400	v3	G	A	70	PASS	.	GT	1/1	0/1	./.
```

## 1) `dnakit.population.fig_voight_painting()`

- **作用：** 围绕单倍型矩阵中的某个变异，绘制 scikit-allel 的双侧 Voight 涂色图，并附带 EHH 衰减曲线。
- **处理方法：** 校验矩阵后调用 `allel.fig_voight_painting()` 并把图保存为 PNG。两侧的涂色数组用 `allel.voight_painting()` 重新计算：左侧作用于 `matrix[:index + 1][::-1]`，右侧作用于 `matrix[index:]`，因此第 0 行是 `index` 处的变异，后面的行逐步远离它。
- **API：** `fig_voight_painting(haplotypes: 'list[list[int]]', output_path: 'str', *, index: 'int | None' = None, palette: 'str' = 'colorblind', height_factor: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`，同时安装 scipy 和 seaborn）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `haplotypes` 是非空的“变异 x 单倍型”矩阵，元素为非负整数（双等位数据用 0 和 1）。`output_path` 必须是新的 `.png` 路径。`index` 是 0 到变异数减 1 之间的整数，默认取中间变异（`n_variants // 2`）。`palette` 是 seaborn 调色板名称，`height_factor` 缩放涂色图高度。还需要 seaborn。
- **输出：** `path`（PNG 路径）和 `data`，其中含 `index`、`n_variants`、`n_haplotypes`、`palette`、`left_painting` 与 `right_painting`（从 `index` 向外的变异 x 单倍型；0 表示该单倍型已不再与其他单倍型共享前缀）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import fig_voight_painting

out = Path(tempfile.mkdtemp())
haplotypes = [
    [0, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1],
]
result = fig_voight_painting(haplotypes, str(out / "voight.png"), index=2, verbose=False)
data = result["data"]
print(data["index"], data["n_variants"], data["n_haplotypes"])
print("left ", data["left_painting"])
print("right", data["right_painting"])
```

- **示例输出：**

```text
2 5 6
left  [[1, 1, 1, 1, 2, 2], [1, 1, 1, 0, 2, 2], [1, 1, 0, 0, 2, 2]]
right [[1, 1, 1, 1, 2, 2], [1, 1, 3, 3, 2, 2], [1, 1, 3, 3, 2, 2]]
```

## 2) `dnakit.population.plot_haplotype_frequencies()`

- **作用：** 把一个区域内各不同单倍型的频率绘制成一条由彩色片段组成的条带。
- **处理方法：** 把矩阵包装成 `allel.HaplotypeArray`，用 `distinct_counts()` 统计相同的单倍型，再用 `allel.plot_haplotype_frequencies()` 绘图。出现多次的单倍型使用调色板颜色，单例使用 `singleton_color`。
- **API：** `plot_haplotype_frequencies(haplotypes: 'list[list[int]]', output_path: 'str', *, palette: 'str' = 'Paired', singleton_color: 'str' = 'w', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`，同时安装 scipy 和 seaborn）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `haplotypes` 是“变异 x 单倍型”的非负整数矩阵；当两列的所有等位基因都相同时它们是同一单倍型。`output_path` 必须是新的 `.png` 路径。`palette`（默认 `Paired`）和 `singleton_color`（默认 `w`，白色）是颜色设置。还需要 seaborn。
- **输出：** `path` 和 `data`，其中含 `n_variants`、`n_haplotypes`、`haplotype_counts`（每种不同单倍型的数量，顺序同 scikit-allel 的 `distinct_counts`）、`singleton_count`、`palette` 和 `singleton_color`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_haplotype_frequencies

out = Path(tempfile.mkdtemp())
haplotypes = [
    [0, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1],
]
result = plot_haplotype_frequencies(haplotypes, str(out / "freq.png"), verbose=False)
data = result["data"]
print(data["haplotype_counts"], data["singleton_count"])
print(data["n_variants"], data["n_haplotypes"])
```

- **示例输出：**

```text
[2, 2, 1, 1] 2
5 6
```

## 3) `dnakit.population.plot_joint_sfs()`

- **作用：** 把两个种群的联合位点频谱绘制成对数刻度图像。
- **处理方法：** 对矩阵调用 `allel.plot_joint_sfs()` 并把图保存为 PNG。默认颜色刻度为对数（`LogNorm`），值为 0 的格子留空。实际绘制的矩阵是输入的转置，所以种群 1 在 x 轴。
- **API：** `plot_joint_sfs(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的二维非负数列表：行是种群 1 的衍生等位基因计数，列是种群 2 的计数。`output_path` 必须是新的 `.png` 路径。`imshow_kwargs` 会传给 matplotlib 的 `imshow`。
- **输出：** `path` 和 `data`，其中含 `spectrum`、`shape`、`drawn_matrix`（实际绘制的转置矩阵），以及坐标轴标签 `xlabel` 与 `ylabel`（各种群的衍生等位基因计数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **示例输出：**

```text
[3, 3] derived allele count (population 1) | derived allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 4) `dnakit.population.plot_joint_sfs_folded()`

- **作用：** 绘制两个种群的折叠联合位点频谱（次要等位基因计数）。
- **处理方法：** 对矩阵调用 `allel.plot_joint_sfs_folded()` 并把图保存为 PNG。颜色刻度与未折叠版本一样是对数的。实际绘制的矩阵是输入的转置，所以种群 1 在 x 轴。
- **API：** `plot_joint_sfs_folded(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的二维非负数列表：行是种群 1 的衍生等位基因计数，列是种群 2 的计数。`output_path` 必须是新的 `.png` 路径。`imshow_kwargs` 会传给 matplotlib 的 `imshow`。
- **输出：** `path` 和 `data`，其中含 `spectrum`、`shape`、`drawn_matrix`（实际绘制的转置矩阵），以及坐标轴标签 `xlabel` 与 `ylabel`（各种群的次要等位基因计数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs_folded

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs_folded(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **示例输出：**

```text
[3, 3] minor allele count (population 1) | minor allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 5) `dnakit.population.plot_joint_sfs_folded_scaled()`

- **作用：** 用 scikit-allel 绘制折叠且缩放的联合位点频谱。
- **处理方法：** 对矩阵调用 `allel.plot_joint_sfs_folded_scaled()` 并把图保存为 PNG。scikit-allel 不会重新缩放数值：按给定值绘制，颜色刻度为线性；请传入 scikit-allel 的 `joint_sfs_folded_scaled` 的结果。实际绘制的矩阵是输入的转置，所以种群 1 在 x 轴。
- **API：** `plot_joint_sfs_folded_scaled(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的二维非负数列表：行是种群 1 的衍生等位基因计数，列是种群 2 的计数。`output_path` 必须是新的 `.png` 路径。`imshow_kwargs` 会传给 matplotlib 的 `imshow`。
- **输出：** `path` 和 `data`，其中含 `spectrum`、`shape`、`drawn_matrix`（实际绘制的转置矩阵），以及坐标轴标签 `xlabel` 与 `ylabel`（各种群的次要等位基因计数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs_folded_scaled

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs_folded_scaled(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **示例输出：**

```text
[3, 3] minor allele count (population 1) | minor allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 6) `dnakit.population.plot_joint_sfs_scaled()`

- **作用：** 用 scikit-allel 的 `plot_joint_sfs_scaled` 绘制缩放后的联合位点频谱。
- **处理方法：** 对矩阵调用 `allel.plot_joint_sfs_scaled()` 并把图保存为 PNG。scikit-allel 不会重新缩放数值：按给定值绘制，颜色刻度为线性，坐标轴标签与未缩放版本相同；请传入 scikit-allel 的 `joint_sfs_scaled` 的结果。实际绘制的矩阵是输入的转置，所以种群 1 在 x 轴。
- **API：** `plot_joint_sfs_scaled(spectrum: 'list[list[float]]', output_path: 'str', *, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的二维非负数列表：行是种群 1 的衍生等位基因计数，列是种群 2 的计数。`output_path` 必须是新的 `.png` 路径。`imshow_kwargs` 会传给 matplotlib 的 `imshow`。
- **输出：** `path` 和 `data`，其中含 `spectrum`、`shape`、`drawn_matrix`（实际绘制的转置矩阵），以及坐标轴标签 `xlabel` 与 `ylabel`（各种群的衍生等位基因计数）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_joint_sfs_scaled

out = Path(tempfile.mkdtemp())
spectrum = [[0, 3, 1], [4, 2, 1], [1, 1, 0]]
result = plot_joint_sfs_scaled(spectrum, str(out / "joint.png"), verbose=False)
data = result["data"]
print(data["shape"], data["xlabel"], "|", data["ylabel"])
print(data["drawn_matrix"])
print(result["path"].endswith("joint.png"))
```

- **示例输出：**

```text
[3, 3] derived allele count (population 1) | derived allele count (population 2)
[[0.0, 4.0, 1.0], [3.0, 2.0, 1.0], [1.0, 1.0, 0.0]]
True
```

## 7) `dnakit.population.plot_moving_haplotype_frequencies()`

- **作用：** 绘制单倍型频率在连续变异窗口中的变化。
- **处理方法：** 把矩阵包装成 `allel.HaplotypeArray` 并调用 `allel.plot_moving_haplotype_frequencies()`。每个窗口的首尾位置由 `allel.moving_statistic()` 计算，排名矩阵则从绘制的网格中读回。图的每一列是含 `window_size` 个变异的窗口，每一行是窗口内单倍型的频率排名（1 为频率最高的非单例，0 为单例）。
- **API：** `plot_moving_haplotype_frequencies(positions: 'list[int]', haplotypes: 'list[list[int]]', window_size: 'int', output_path: 'str', *, start: 'int' = 0, stop: 'int | None' = None, n: 'int | None' = None, palette: 'str' = 'Paired', singleton_color: 'str' = 'w', verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`，同时安装 scipy 和 seaborn）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `positions` 是非负整数列表，与 `haplotypes` 的行一一对应。`haplotypes` 是“变异 x 单倍型”的非负整数矩阵。`window_size` 是不小于 1 的整数。`output_path` 必须是新的 `.png` 路径。`start` 与 `stop` 是变异索引，`n`（不小于 1）只保留频率最高的 n 个排名。还需要 seaborn。
- **输出：** `path` 和 `data`，其中含 `window_size`、`start`、`stop`、`n`、`window_bounds`（每个窗口的首尾位置）、`rank_matrix`（窗口 x 单倍型，即实际绘制的矩阵）、`palette` 和 `singleton_color`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_moving_haplotype_frequencies

out = Path(tempfile.mkdtemp())
positions = [10, 25, 40, 70, 90, 120, 150, 170]
haplotypes = [
    [0, 0, 1, 1, 0, 1],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 1, 1, 0, 0],
    [0, 0, 1, 1, 1, 1],
]
result = plot_moving_haplotype_frequencies(
    positions, haplotypes, 4, str(out / "moving.png"), verbose=False
)
data = result["data"]
print(data["window_bounds"])
print(data["rank_matrix"])
top = plot_moving_haplotype_frequencies(
    positions, haplotypes, 4, str(out / "top1.png"), n=1, verbose=False
)
print(top["data"]["rank_matrix"])
```

- **示例输出：**

```text
[[10, 70], [90, 170]]
[[1, 1, 2, 2, 0, 0], [1, 1, 2, 2, 3, 3]]
[[1, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0]]
```

## 8) `dnakit.population.plot_pairwise_distance()`

- **作用：** 把样本间的两两距离矩阵绘制成图像，可带样本标签。
- **处理方法：** 检查输入是方阵或含 n*(n-1)/2 个值的压缩向量，然后调用 `allel.plot_pairwise_distance()`，它会把压缩向量展开为对角线为零的方阵。矩阵从绘制的图像中读回。
- **API：** `plot_pairwise_distance(distances: 'list[float]', output_path: 'str', *, labels: 'list[str] | None' = None, colorbar: 'bool' = True, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `distances` 是有限数值列表：要么是按行展开的上三角压缩形式，要么是 n x n 的方阵（其他长度会引发 `ConfigurationError`）。`labels` 是 n 个字符串的列表。`colorbar` 控制是否添加颜色条。`output_path` 必须是新的 `.png` 路径。
- **输出：** `path` 和 `data`，其中含 `n_samples`、`labels`、`colorbar` 和 `square_matrix`（实际绘制的矩阵）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_pairwise_distance

out = Path(tempfile.mkdtemp())
condensed = [0.1, 0.2, 0.3]  # (a,b), (a,c), (b,c)
result = plot_pairwise_distance(
    condensed, str(out / "distance.png"), labels=["a", "b", "c"], verbose=False
)
data = result["data"]
print(data["n_samples"], data["labels"])
print(data["square_matrix"])
```

- **示例输出：**

```text
3 ['a', 'b', 'c']
[[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]]
```

## 9) `dnakit.population.plot_pairwise_ld()`

- **作用：** 把所有变异对之间的连锁不平衡（r 平方）绘制成三角形图像。
- **处理方法：** 检查输入是方阵或含 n*(n-1)/2 个值的压缩向量，然后调用 `allel.plot_pairwise_ld()`。scikit-allel 保留方阵的下三角并上下翻转；矩阵从绘制的图像中读回。该函数只绘制传入的数值，不计算 LD。
- **API：** `plot_pairwise_ld(ld_values: 'list[float]', output_path: 'str', *, colorbar: 'bool' = True, imshow_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `ld_values` 是有限数值列表：按行展开的上三角压缩形式或方阵。`colorbar` 控制是否添加颜色条，`imshow_kwargs` 传给 matplotlib 的 `imshow`。`output_path` 必须是新的 `.png` 路径。
- **输出：** `path` 和 `data`，其中含 `n_variants`、`ld_values`（输入）、`colorbar` 和 `drawn_matrix`（实际绘制的翻转后的下三角）。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_pairwise_ld

out = Path(tempfile.mkdtemp())
condensed = [0.5, 0.1, 0.9]  # r squared of pairs (0,1), (0,2), (1,2)
result = plot_pairwise_ld(condensed, str(out / "ld.png"), verbose=False)
data = result["data"]
print(data["n_variants"])
print(data["drawn_matrix"])
```

- **示例输出：**

```text
3
[[0.1, 0.9, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]]
```

## 10) `dnakit.population.plot_sfs()`

- **作用：** 绘制一维位点频谱。
- **处理方法：** 对频谱调用 `allel.plot_sfs()` 并把图保存为 PNG；点坐标从绘制的曲线读回。
- **API：** `plot_sfs(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的一维非负位点计数列表，按等位基因计数索引。`output_path` 必须是新的 `.png` 路径。`yscale` 是 matplotlib 刻度名（默认 `log`）。`bins` 是整数或分箱边界列表。`n` 是采样的染色体数，设置后 x 轴变为等位基因频率。`clip_endpoints`（默认 true）去掉首尾两项。`label` 与 `plot_kwargs` 传给曲线。
- **输出：** `path` 和 `data`，其中含 `spectrum`、绘制的点 `x` 与 `y`、`xlabel`（`derived allele count`；给定 `n` 时为 `derived allele frequency`）、`ylabel`（`site frequency`）、`yscale`、`n`、`bins`、`clip_endpoints` 和 `label`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
with_n = plot_sfs(spectrum, str(out / "sfs_n.png"), n=6, yscale="linear", verbose=False)
print(with_n["data"]["x"][:2], with_n["data"]["xlabel"])
```

- **示例输出：**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
derived allele count | site frequency | log
[0.16666666666666666, 0.3333333333333333] derived allele frequency
```

## 11) `dnakit.population.plot_sfs_folded()`

- **作用：** 绘制折叠的一维位点频谱（次要等位基因计数）。
- **处理方法：** 对频谱调用 `allel.plot_sfs_folded()` 并把图保存为 PNG；点坐标从绘制的曲线读回。参数与 `plot_sfs` 相同，只有 x 轴标签不同。
- **API：** `plot_sfs_folded(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的一维非负位点计数列表，按等位基因计数索引。`output_path` 必须是新的 `.png` 路径。`yscale` 是 matplotlib 刻度名（默认 `log`）。`bins` 是整数或分箱边界列表。`n` 是采样的染色体数，设置后 x 轴变为等位基因频率。`clip_endpoints`（默认 true）去掉首尾两项。`label` 与 `plot_kwargs` 传给曲线。
- **输出：** `path` 和 `data`，其中含 `spectrum`、绘制的点 `x` 与 `y`、`xlabel`（`minor allele count`）、`ylabel`（`site frequency`）、`yscale`、`n`、`bins`、`clip_endpoints` 和 `label`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_folded

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs_folded(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
```

- **示例输出：**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
minor allele count | site frequency | log
```

## 12) `dnakit.population.plot_sfs_folded_scaled()`

- **作用：** 绘制折叠且缩放的一维位点频谱。
- **处理方法：** 对频谱调用 `allel.plot_sfs_folded_scaled()` 并把图保存为 PNG；点坐标从绘制的曲线读回。scikit-allel 按给定值绘制，仅改变标签；要绘制缩放后的值，请传入其 `sfs_folded_scaled` 的结果。
- **API：** `plot_sfs_folded_scaled(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的一维非负位点计数列表，按等位基因计数索引。`output_path` 必须是新的 `.png` 路径。`yscale` 是 matplotlib 刻度名（默认 `log`）。`bins` 是整数或分箱边界列表。`n` 是采样的染色体数，设置后 x 轴变为等位基因频率。`clip_endpoints`（默认 true）去掉首尾两项。`label` 与 `plot_kwargs` 传给曲线。
- **输出：** `path` 和 `data`，其中含 `spectrum`、绘制的点 `x` 与 `y`、`xlabel`（`minor allele count`）、`ylabel`（`scaled site frequency`）、`yscale`、`n`、`bins`、`clip_endpoints` 和 `label`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_folded_scaled

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs_folded_scaled(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
```

- **示例输出：**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
minor allele count | scaled site frequency | log
```

## 13) `dnakit.population.plot_sfs_scaled()`

- **作用：** 绘制缩放后的一维位点频谱。
- **处理方法：** 对频谱调用 `allel.plot_sfs_scaled()` 并把图保存为 PNG；点坐标从绘制的曲线读回。scikit-allel 按给定值绘制，仅改变 y 轴标签；要绘制缩放后的值，请传入其 `sfs_scaled` 的结果。
- **API：** `plot_sfs_scaled(spectrum: 'list[float]', output_path: 'str', *, yscale: 'str' = 'log', bins: 'int | list[float] | None' = None, n: 'int | None' = None, clip_endpoints: 'bool' = True, label: 'str | None' = None, plot_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `spectrum` 是非空的一维非负位点计数列表，按等位基因计数索引。`output_path` 必须是新的 `.png` 路径。`yscale` 是 matplotlib 刻度名（默认 `log`）。`bins` 是整数或分箱边界列表。`n` 是采样的染色体数，设置后 x 轴变为等位基因频率。`clip_endpoints`（默认 true）去掉首尾两项。`label` 与 `plot_kwargs` 传给曲线。
- **输出：** `path` 和 `data`，其中含 `spectrum`、绘制的点 `x` 与 `y`、`xlabel`（`derived allele count`）、`ylabel`（`scaled site frequency`）、`yscale`、`n`、`bins`、`clip_endpoints` 和 `label`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_scaled

out = Path(tempfile.mkdtemp())
spectrum = [0, 40, 20, 10, 5, 2, 0]
result = plot_sfs_scaled(spectrum, str(out / "sfs.png"), verbose=False)
data = result["data"]
print(data["x"], data["y"])
print(data["xlabel"], "|", data["ylabel"], "|", data["yscale"])
```

- **示例输出：**

```text
[1, 2, 3, 4, 5] [40.0, 20.0, 10.0, 5.0, 2.0]
derived allele count | scaled site frequency | log
```

## 14) `dnakit.population.plot_variant_locator()`

- **作用：** 绘制从每个变异的基因组位置到其在变异索引中位置的连线。
- **处理方法：** 检查位置已排序后调用 `allel.plot_variant_locator()`。未给出 `step` 时，DNAKit 使用 `max(1, n_variants // 100)`（变异少于 100 个时 scikit-allel 自身会因步长为 0 而失败）。x 轴范围和绘制的线数从坐标轴读取。
- **API：** `plot_variant_locator(positions: 'list[int]', output_path: 'str', *, step: 'int | None' = None, start: 'int | None' = None, stop: 'int | None' = None, flip: 'bool' = False, line_kwargs: 'dict[str, Any] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `positions` 是单个 contig 的非负整数列表，必须升序（未排序会引发 `ConfigurationError`）。`step` 是不小于 1 的整数。`start` 与 `stop` 限定基因组区域。`flip` 翻转图形，`line_kwargs` 传给连线。`output_path` 必须是新的 `.png` 路径。
- **输出：** `path` 和 `data`，其中含 `n_variants`、`step`、`start` 与 `stop`（实际绘制的 x 轴范围）、`flip` 和 `lines_drawn`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_variant_locator

out = Path(tempfile.mkdtemp())
positions = [10, 25, 40, 70, 90, 120, 150, 170]
result = plot_variant_locator(
    positions, str(out / "locator.png"), step=2, start=0, stop=200, verbose=False
)
data = result["data"]
print(data["n_variants"], data["step"], data["start"], data["stop"], data["lines_drawn"])
default = plot_variant_locator(positions, str(out / "locator2.png"), verbose=False)
print(default["data"]["step"], default["data"]["lines_drawn"])
```

- **示例输出：**

```text
8 2 0.0 200.0 4
1 8
```

## 15) `dnakit.population.plot_voight_painting()`

- **作用：** 绘制单侧 Voight 涂色图，其中每个共享的单倍型前缀是一块彩色区域。
- **处理方法：** 校验涂色数组后调用 `allel.plot_voight_painting()`；`flank='left'` 会反转变异轴。颜色数是数组中最大的分组编号。
- **API：** `plot_voight_painting(painting: 'list[list[int]]', output_path: 'str', *, palette: 'str' = 'colorblind', flank: 'str' = 'right', height_factor: 'float' = 0.01, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`，同时安装 scipy 和 seaborn）和 matplotlib（extra `dnakit[bio-plots]`）；离线运行。
- **输入：** `painting` 是非空的二维非负整数列表，通常是 scikit-allel 的 `voight_painting` 的第一个返回值（变异 x 单倍型；0 表示单例）。`flank` 为 `right` 或 `left`，其他值会引发 `ConfigurationError`。`palette` 与 `height_factor` 同 `fig_voight_painting`。`output_path` 必须是新的 `.png` 路径。还需要 seaborn。
- **输出：** `path` 和 `data`，其中含 `painting`、`flank`、`n_colors`、`palette` 和 `height_factor`。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_voight_painting

out = Path(tempfile.mkdtemp())
painting = [
    [1, 1, 2, 2, 0],
    [1, 1, 2, 3, 0],
    [1, 1, 0, 0, 0],
]
result = plot_voight_painting(painting, str(out / "painting.png"), flank="left", verbose=False)
data = result["data"]
print(data["flank"], data["n_colors"])
try:
    plot_voight_painting(painting, str(out / "bad.png"), flank="up", verbose=False)
except Exception as error:
    print(type(error).__name__, error.message)
```

- **示例输出：**

```text
left 3
ConfigurationError flank must be 'right' or 'left'.
```

## 16) `dnakit.population.roh_poissonhmm()`

- **作用：** 用 scikit-allel 的泊松隐马尔可夫模型为单个个体调用纯合片段（ROH）。
- **处理方法：** 校验基因型和位置、检查窗口与 contig 设置，要求 `protopunica` 包可用，然后用 ROH 内的杂合率（`phet_roh`）和 ROH 外的两个杂合率（`phet_nonroh_low`、`phet_nonroh_high`）调用 `allel.roh_poissonhmm()`。ROH 表以记录列表返回，并给出 FROH，即 ROH 覆盖 contig 的比例。
- **API：** `roh_poissonhmm(genotypes: 'list[list[int]]', positions: 'list[int]', *, phet_roh: 'float' = 0.001, phet_nonroh_low: 'float' = 0.0025, phet_nonroh_high: 'float' = 0.01, transition: 'float' = 0.001, window_size: 'int' = 1000, min_roh: 'int' = 0, is_accessible: 'list[bool] | None' = None, contig_size: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）以及 `protopunica`。`protopunica` 不在任何 DNAKit extra 中，需要单独用 pip 安装；缺少时抛出带安装提示的 `BackendUnavailableError`。本环境未安装 `protopunica`，所以下面的示例展示参数检查和缺少依赖时的真实行为。
- **输入：** `genotypes` 是“变异 x 倍性”的非负等位基因索引列表（二倍体为两列），`positions` 是对应的 1 起始位置。必须给出 `is_accessible`（contig 每个碱基一个真/假值）或 `contig_size`，否则引发 `ConfigurationError`。`window_size` 不小于 1，`min_roh` 不小于 0。
- **输出：** 字典，含 `status`、`roh_blocks`（来自 scikit-allel ROH 表的记录列表）、`block_count`、`froh`、`window_size`、`phet_roh`、`phet_nonroh`、`transition` 和 `min_roh`。
- **示例代码：**

```python
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import roh_poissonhmm

genotypes = [[0, 0], [0, 1], [1, 1]]
positions = [100, 200, 300]
try:
    roh_poissonhmm(genotypes, positions, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, error.message)
try:
    roh_poissonhmm(genotypes, positions, contig_size=1000, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, "protopunica could not be imported" in error.message)
```

- **示例输出：**

```text
ConfigurationError Provide is_accessible (one flag per base) or contig_size.
BackendUnavailableError True
```

安装依赖后的用法（Usage with the dependency installed，不执行）：

```python
from dnakit.population import roh_poissonhmm

# genotypes: one row per variant, two allele columns for a diploid individual
result = roh_poissonhmm(
    genotypes,
    positions,
    contig_size=5_000_000,
    window_size=1000,
    verbose=False,
)
print(result["block_count"], result["froh"])
```

## 17) `dnakit.population.allel_vcf_to_zarr()`

- **作用：** 用 scikit-allel 读取 VCF 文件，并把变异与基因型数组写入新的 Zarr 存储。
- **处理方法：** 用 `allel.read_vcf()` 读取文件，并用 zarr 3 的 `create_array()` 按 `vcf_to_zarr` 的布局（`variants/<FIELD>`、`calldata/<FIELD>`、`samples`）写出每个数组。它不调用 scikit-allel 自带的 `vcf_to_zarr`，因为后者使用 zarr 3 已没有的 zarr 2 接口。字符串数组按字符串存储。VCF 没有变异时不写存储。
- **API：** `allel_vcf_to_zarr(vcf_path: 'str | Path', zarr_path: 'str | Path', *, fields: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`。
- **运行要求：** 需要可选包 scikit-allel（extra `dnakit[bio-scikit-allel]`）和 zarr 3（extra `dnakit[io-zarr]`）；离线运行。注意：缺少 zarr 时错误提示给出的 extra 是 `bio-scikit-allel`，而该 extra 并不包含 zarr。
- **输入：** `vcf_path` 必须是已存在的文件。`zarr_path` 必须尚不存在，且其父目录必须存在，否则引发 `ConfigurationError`。`fields` 是最多 200 个 scikit-allel 字段名（如 `variants/CHROM`、`calldata/GT`）的非空列表，或 `None` 表示使用 scikit-allel 的默认字段。整个 VCF 会被读入内存。
- **输出：** 字典，含 `status`、`operation`、`output_path`（存储的绝对路径；无变异的 VCF 为 `None`）、`variants`、`samples`、`arrays`（每个数组的名称、形状与 dtype）和 `provenance`（包版本与布局）。无变异的 VCF 返回 `variants` 为 0、`arrays` 为空列表，且没有 `samples` 键。
- **示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.population import allel_vcf_to_zarr

out = Path(tempfile.mkdtemp())
result = allel_vcf_to_zarr("mini.vcf", out / "mini.zarr", verbose=False)
print(result["status"], result["variants"], result["samples"])
for array in result["arrays"]:
    print(array["name"], array["shape"])
only_gt = allel_vcf_to_zarr(
    "mini.vcf", out / "gt.zarr", fields=["calldata/GT"], verbose=False
)
print([array["name"] for array in only_gt["arrays"]])
try:
    allel_vcf_to_zarr("mini.vcf", out / "mini.zarr", verbose=False)
except Exception as error:
    print(type(error).__name__, "already exists" in error.message)
```

- **示例输出：**

```text
completed 3 3
calldata/GT [3, 3, 2]
samples [3]
variants/ALT [3, 3]
variants/CHROM [3]
variants/FILTER_PASS [3]
variants/ID [3]
variants/POS [3]
variants/QUAL [3]
variants/REF [3]
['calldata/GT']
ConfigurationError True
```
