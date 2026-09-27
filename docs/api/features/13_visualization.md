# 可视化

绘制序列、注释、比对和样本分布。PNG/JPG 导出需要 `viz` 依赖，见 [FAQ：补充说明](../../faq/tutorial_basics.md#visualization-export)。

## 1) 序列文字图 {#1-viz-001}

- **作用：** 把 DNA 碱基、坐标和可选注释绘制成矢量 SVG，返回可嵌入网页或保存的图形，用于快速查看短序列内容。
- **计算方法：** 将碱基按每行数量排版，按位置绘制字符、坐标与可选互补链。
- **API：** `dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`；`config` 使用 `dnakit.visualization.SequencePlotConfig`。
- **输入：** `DNASequence` 或 `DNARecord`；可选每行碱基数、坐标、互补链、行间距、列间距、自定义字符映射和显示上限。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.visualization import SequencePlotConfig, plot_sequence, save_image

artifact = plot_sequence(
    DNASequence("ACGTGGTCCA"),
    config=SequencePlotConfig(bases_per_line=4),
)
save_image(artifact, "viz-001-sequence", image_type="png")
```

- **示例结果（`viz-001-sequence.png`）：**

![序列文字图示例](../../assets/images/visualization/viz-001-sequence.png){ width="420" }

## 2) 位置高亮 {#2-viz-002}

- **作用：** 按坐标、颜色和标签高亮指定碱基或区间，使 motif、突变和其他关注位置在序列图中清晰可见。
- **计算方法：** 把半开区间映射到字符位置，再叠加颜色、标签与优先级。
- **API：** `dnakit.visualization.Highlight(start[必须], end[必须], label[可选], color[可选], foreground[可选], priority[可选], opacity[可选])`、`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`。
- **输入：** 序列和高亮区间；可选标签、颜色和优先级。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.visualization import Highlight, plot_sequence, save_image

artifact = plot_sequence(
    DNASequence("ACGTACGT"),
    highlights=[Highlight(1, 4, label="target", color="#ffcc00")],
)
save_image(artifact, "viz-002-highlight", image_type="png")
```

- **示例结果（`viz-002-highlight.png`）：**

![位置高亮示例](../../assets/images/visualization/viz-002-highlight.png){ width="420" }

## 3) 样式控制 {#3-viz-003}

- **作用：** 通过统一样式配置控制颜色、字体、尺寸、边距、行间距、文字列间距和布局，使不同可视化结果具有一致外观并可复现。
- **计算方法：** 按 SequencePlotConfig 的字体、间距、列数和颜色生成 SVG 元素。
- **API：** `dnakit.visualization.SVGTheme(background[可选], foreground[可选], muted[可选], grid[可选], accent[可选], missing[可选], font_family[可选])`；序列文字图使用 `SequencePlotConfig(font_size[可选], cell_width[可选], line_height[可选], column_spacing[可选], line_spacing[可选], margin[可选], theme[可选])`。
- **输入：** 绘图配置；`column_spacing` 控制相邻文字列之间增加的水平像素，`line_spacing` 控制相邻序列行之间增加的垂直像素，二者默认均为 `0`。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.visualization import SVGTheme, SequencePlotConfig, plot_sequence, save_image

theme = SVGTheme(background="#111827", foreground="#f9fafb", accent="#22c55e")
artifact = plot_sequence(
    DNASequence("ACGTGGTCCA"),
    config=SequencePlotConfig(
        bases_per_line=4,
        font_size=18,
        column_spacing=8,
        line_spacing=18,
        theme=theme,
    ),
)
save_image(artifact, "viz-003-theme", image_type="png")
```

- **示例结果（`viz-003-theme.png`）：**

![深色样式控制示例](../../assets/images/visualization/viz-003-theme.png){ width="420" }

## 4) Gap 显示 {#4-viz-004-gap}

- **作用：** 在序列图中用专门符号显示 Gap 的位置、已知或未知长度，同时保持坐标含义，避免把缺失区域误画成真实碱基。
- **计算方法：** 根据 Gap 长度信息绘制缺口标记，保留已知与未知缺口的区别。
- **API：** `dnakit.Gap(length[必须], kind[可选], crossable[可选], evidence[可选], metadata[可选])`、`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`。
- **输入：** 包含显式 `Gap` part 的 `DNASequence` 或记录。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNASequence, Gap
from dnakit.visualization import plot_sequence, save_image

sequence = DNASequence(["AC", Gap(500), "T", Gap(None), "G"])
artifact = plot_sequence(sequence)
save_image(artifact, "viz-004-gap", image_type="png")
```

- **示例结果（`viz-004-gap.png`）：**

![Gap 显示示例](../../assets/images/visualization/viz-004-gap.png){ width="432" }

## 5) 线性序列图 {#5-viz-005}

- **作用：** 沿线性 DNA 坐标绘制 feature 的起止范围、方向和标签，生成类似基因结构图的 SVG，用于查看多个注释的相对位置。
- **计算方法：** 把 feature 区间线性映射到画布坐标，再按方向绘制注释。
- **API：** `dnakit.visualization.plot_linear_map(value[必须], width[可选], height[可选], max_features[可选], title[可选], theme[可选])`；画布边长取 `width` 和 `height` 中较大的值，`title` 仅作为 SVG 无障碍说明。
- **输入：** resolved linear `DNASequence` 或带 `DNAFeature` 的 `DNARecord`；输出始终为正方形。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNAFeature, DNARecord, DNASequence, Interval
from dnakit.visualization import plot_linear_map, save_image

record = DNARecord(
    DNASequence("ACGTACGT"),
    "linear",
    features=[DNAFeature("gene", Interval(1, 6))],
)
artifact = plot_linear_map(record)
save_image(artifact, "viz-005-linear-map", image_type="png")
```

- **示例结果（`viz-005-linear-map.png`）：**

![线性序列图示例](../../assets/images/visualization/viz-005-linear-map.png){ width="900" }

## 6) 环状 DNA 图 {#6-viz-006-dna}

- **作用：** 把环状 DNA 及其 feature 映射到圆周坐标，显示方向、跨原点区间和标签，用于质粒或其他环状分子的结构概览。
- **计算方法：** 将序列位置转换为圆周角度，将 feature 区间绘制为定向圆弧。
- **API：** `dnakit.visualization.plot_circular_map(value[必须], size[可选], max_features[可选], title[可选], theme[可选])`；`size` 是正方形画布边长，`title` 仅作为 SVG 无障碍说明。
- **输入：** `topology="circular"` 的 resolved 序列或记录。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNAFeature, DNARecord, DNASequence, Interval
from dnakit.visualization import plot_circular_map, save_image

record = DNARecord(
    DNASequence("ACGTACGT", topology="circular"),
    "plasmid",
    features=[DNAFeature("gene", Interval(1, 6))],
)
artifact = plot_circular_map(record)
save_image(artifact, "viz-006-circular-map", image_type="png")
```

- **示例结果（`viz-006-circular-map.png`）：**

![环状 DNA 图示例](../../assets/images/visualization/viz-006-circular-map.png){ width="520" }

## 7) 自定义字符表示 {#7-viz-007}

- **作用：** 把 A、T、C、G 或 IUPAC 符号替换为调用方指定的显示字符，例如 A 显示为 `*`、T 显示为 `-`、C 显示为 `+`、G 显示为 `]`；坐标、高亮、互补关系和碱基颜色仍按原始序列计算。
- **计算方法：** 先按原始碱基计算坐标与配色，再用 symbol_map 替换显示字符。
- **API：** `dnakit.visualization.SequencePlotConfig(symbol_map[可选])`；通过 `plot_sequence(value, config=config)` 使用。
- **输入：** `symbol_map` 的键是大写 DNA/IUPAC 符号，值是一个可见字符；未配置的符号保持原样，默认映射为空。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.visualization import SequencePlotConfig, plot_sequence, save_image

artifact = plot_sequence(
    DNASequence("ATCGATCG"),
    config=SequencePlotConfig(
        bases_per_line=4,
        font_size=20,
        column_spacing=8,
        line_spacing=14,
        symbol_map={"A": "*", "T": "-", "C": "+", "G": "]"},
    ),
)
save_image(artifact, "viz-007-custom-symbols", image_type="png")
```

- **示例结果（`viz-007-custom-symbols.png`）：**

![自定义字符表示示例](../../assets/images/visualization/viz-007-custom-symbols.png){ width="420" }

## 8) Alignment 图 {#8-viz-008-alignment}

- **作用：** 将成对比对的两条序列按列绘制，区分匹配、错配、插入和删除，并显示坐标，便于人工检查差异位置。
- **计算方法：** 按比对结果逐列绘制两个序列和 Gap，并标识匹配与错配。
- **API：** `dnakit.visualization.plot_alignment(result[必须], columns_per_line[可选], max_columns[可选], theme[可选])`。
- **输入：** 预先计算的 `AlignmentResult`；可选每行列数和上限，输出始终为正方形。
- **输出：** 返回 `SVGArtifact` 图形；示例中的 `save_image` 另外将图形保存为图片文件。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.alignment import align_pairwise
from dnakit.visualization import plot_alignment, save_image

alignment = align_pairwise(DNASequence("ACGT"), DNASequence("AGT"))
artifact = plot_alignment(alignment, columns_per_line=4)
save_image(artifact, "viz-008-alignment", image_type="png")
```

- **示例结果（`viz-008-alignment.png`）：**

![Alignment 图示例](../../assets/images/visualization/viz-008-alignment.png){ width="640" }

## 9) 分类与回归样本分布 {#viz-009}

- **作用：** 按类别或数值区间绘制样本占比图，也可直接读取 CSV 指定列。
- **计算方法：** 每根柱的高度为该类别或区间的有效样本数除以该列有效样本总数；所有占比之和为 1。
- **API：** `dnakit.visualization.plot_distribution(values, task=..., config=...)`；读取 CSV 使用 `plot_csv_distributions(source, columns=..., task=..., config=...)`。
- **真实与预测对比：** `plot_label_comparison(true_labels, predicted_labels, config=..., output=..., overwrite=...)` 接收两个等长的标签列表，按类别并排绘制真实与预测的占比。`DistributionPlotConfig(categories=("0", "1"), show_values=True)` 可固定类别顺序并显示柱顶百分比；传入 `output="comparison.svg"` 即可在绘图时保存，格式由文件后缀决定，默认不覆盖已有文件。返回的 `metadata` 分别包含两组计数与占比。空值分别从各组的分母中排除。
- **输入：** 类别标签、连续数值或带表头的 CSV。`probability` 模式要求值在 `[0,1]`；数值模式拒绝 NaN/Inf。
- **输出：** 返回 `SVGArtifact`；`metadata` 中包含各类别或区间的数量与占比。
- **示例代码：**

```python
from dnakit.visualization import plot_distribution

result = plot_distribution([0, 1, 1, 2, 2], task="classification")
print(result.metadata["counts"])
print(result.metadata["proportions"])

from dnakit.visualization import DistributionPlotConfig, plot_label_comparison

comparison = plot_label_comparison(
    [0, 0, 1, 1], [0, 1, 1, 1],
    config=DistributionPlotConfig(categories=("0", "1"), show_values=True),
    output="label_comparison.svg",
    overwrite=True,
)
print(comparison.metadata["true_proportions"])
print(comparison.metadata["predicted_proportions"])
```

- **示例结果：**

```text
(1, 2, 2)
(0.2, 0.4, 0.4)
(0.5, 0.5)
(0.25, 0.75)
```

[FAQ：补充说明](../../faq/tutorial_basics.md#distribution-plots)

### 从 CSV 绘制指标对比图

`plot_metric_comparison_csv` 读取包含 `model,mean` 的 CSV；可选列 `se` 和 `source` 分别提供标准误与分组。只需 CSV 路径和图片输出路径即可绘图，文件后缀决定 SVG、PNG 或 PDF 格式。绘图需要安装 `bio-plots` 可选依赖。还可通过 `fold_csvs` 读取逐折指标 CSV 并计算均值及标准误，通过 `test_csv` 加入测试集指标；不同来源显示在独立面板。

```csv
source,model,mean,se
Reference,Model A,0.614,0.010
Reference,Model B,0.517,0.011
```

```python
from dnakit.visualization import plot_metric_comparison_csv

plot_metric_comparison_csv("metrics.csv", output="comparison.svg", metric_label="MCC")
```

### 从预测 CSV 绘制混淆矩阵

`plot_confusion_matrix_csv` 只需一份包含真实标签 `label` 和预测标签 `prediction` 的 CSV。行表示真实类别，列表示预测类别；左图显示样本数，右图显示每个真实类别中各预测类别的比例。`labels` 可固定类别顺序，也可保留本次数据未出现的类别。传入 `output` 会直接保存图片；SVG 无需额外绘图依赖，PNG/PDF 导出需要 `viz` 可选依赖。

```csv
label,prediction
0,0
0,1
1,1
```

```python
from dnakit.visualization import plot_confusion_matrix_csv

plot_confusion_matrix_csv(
    "predictions.csv",
    labels=("0", "1"),
    output="confusion.svg",
    overwrite=True,
)
```

如果 CSV 使用不同的列名，可通过 `true_column` 和 `predicted_column` 指定。函数返回的 `SVGArtifact.metadata` 包含 `counts`、`row_normalized` 和 `row_totals`。

### 从逐轮 CSV 绘制学习曲线

`plot_learning_curves_csv` 从一份 CSV 绘制训练集和验证集的 loss、MCC 曲线，并直接保存图片。必需列为 `epoch,train_loss,validation_loss,train_mcc,validation_mcc`；可选的 `dataset,setting,method,fold` 分别表示数据集、配置、方法和交叉验证折。相同方法在同一轮的多折结果取均值，阴影表示样本标准差；空的指标单元格会跳过。若 CSV 包含多个数据集或配置，传入输出目录即可为每组生成一张 SVG；仅有一组时也可直接传入 `.svg`、`.png` 或 `.pdf` 图片路径。此功能需要 `bio-plots` 可选依赖。

```csv
dataset,setting,method,fold,epoch,train_loss,validation_loss,train_mcc,validation_mcc
enhancers,baseline,full,0,1,0.54,0.50,0.59,0.52
enhancers,baseline,full,1,1,0.56,0.48,0.57,0.54
```

```python
from dnakit.visualization import plot_learning_curves_csv

plot_learning_curves_csv("epochs.csv", output="figures", overwrite=True)
```

上述示例生成 `figures/enhancers_baseline.svg`；不提供分组列时，一份 CSV 生成一张图。

### 从 DNA embedding 和错误标记 CSV 绘制误判投影

`plot_embedding_errors` 读取一个 NPZ 表征文件和一个带 `id,error` 列的 CSV，按样本 ID 对齐后将表征投影到二维，并高亮错误预测。NPZ 必须包含 `embeddings`（二维数值矩阵）和 `ids`；可选 `splits`（`train`/`test`）和 `labels`。有训练集时，仅用训练集拟合标准化和 PCA，绘制测试集；没有划分时绘制全部样本。CSV 必须覆盖每个要绘制的 ID，`error` 可为 `correct`、`FP`、`FN`、`error` 或 `0`/`1`；可选 `method` 列可在一个 CSV 中表示多种模型，每种模型生成一张 SVG。`label`、`prediction`、`probability` 列可选；提供预测列时会核对错误标记。需要 `bio-analysis` 和 `bio-plots` 依赖。

```csv
id,method,error,label,prediction
sample-1,full,correct,0,0
sample-2,full,FN,1,0
```

```python
from dnakit.visualization import plot_embedding_errors

plots = plot_embedding_errors("embeddings.npz", "predictions.csv")
# 默认写入 predictions/full.svg；可用 output_dir 指定目录。
```

可选 `n_clusters=5, save_details=True` 会保存 PCA 投影、逐点 CSV 和各簇错误率 CSV。`overwrite=True` 可覆盖已有输出。

<span id="_2"></span>
<span id="_3"></span>
<span id="csv"></span>
<span id="python-api"></span>

## 10) 绘制进化树 {#bio-7-plot}

- **作用：** 查看树的枝长和叶节点，并保存为图片。
- **计算方法：** 使用 Matplotlib 绘制树的分支与标签。
- **API：** `dnakit.bioinformatics.phylogeny.plot_tree(tree)`
- **输入：** `TreeNode` 树对象；需要 `bio-plots`。示例自行创建树。
- **输出：** 返回 Matplotlib `Axes`；可通过其 `figure` 保存树图。
- **示例代码：**

```python
from dnakit.bioinformatics.phylogeny import TreeNode, plot_tree

tree = TreeNode(children=(TreeNode("a", 1), TreeNode("b", 1)))
axes = plot_tree(tree)
print(type(axes).__name__)
axes.figure.savefig("tree.png")
```

- **示例结果：**

```text
Axes
```

树图保存为 `tree.png`。

[依赖与实现说明](../../faq/bioinformatics.md#alignment)。

## 11) 绘制序列 dotplot {#bio-51}

- **作用：** 查看两条序列或单条序列内部的重复和反向互补匹配。
- **计算方法：** 查找正反链精确 k-mer 匹配并绘制坐标点。
- **API：** `dnakit.bioinformatics.motifs.dotplot(left, right=None, word_length=10, max_matches=1000000)`
- **输入：** 一条或两条序列、匹配词长；需要 `bio-plots`。
- **输出：** 返回 Matplotlib `Axes`，展示序列匹配点图。
- **示例代码：**

```python
from dnakit.bioinformatics.motifs import dotplot

axes = dotplot("ACGTACGT", "ACGTTCGT", word_length=3)
print(type(axes).__name__)
# axes.figure.savefig("dotplot.png")
```

- **示例结果：**

```text
Axes
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 12) 绘制 Sanger 色谱 {#bio-55-plot}

- **作用：** 查看 Sanger 测序各碱基通道的信号，并保存色谱图。
- **计算方法：** 使用 Matplotlib 绘制 ABI 信号通道及碱基标注。
- **API：** `dnakit.bioinformatics.cloning.plot_chromatogram(trace)`
- **输入：** `read_sanger` 返回的 trace 字典；示例读取已有 `read.ab1`，需要 `bio-formats` 和 `bio-plots`。
- **输出：** 返回 Matplotlib `Axes`，可保存为 Sanger 色谱图。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import read_sanger, plot_chromatogram

trace = read_sanger("read.ab1")
axes = plot_chromatogram(trace)
axes.figure.savefig("chromatogram.png")
```

- **示例结果：**

色谱图保存为 `chromatogram.png`。需要真实 ABI，此处未执行。

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 13) 绘制定性凝胶条带图 {#bio-58}

- **作用：** 根据片段长度预览多泳道的条带分布。
- **计算方法：** 使用 log(片段长度) 对应迁移位置；各片段按相同上样量绘制。
- **API：** `dnakit.bioinformatics.cloning.simulate_gel(lanes, min_bp=100, max_bp=10000)`
- **输入：** 每泳道的长度列表或线性 `Duplex` 列表；需要 `bio-plots`。用于定性示意。
- **输出：** 返回 Matplotlib `Axes`，展示各泳道的模拟条带。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import simulate_gel

axes = simulate_gel([[500, 1000, 2000], [750, 1500]])
print(type(axes).__name__)
# axes.figure.savefig("gel.png")
```

- **示例结果：**

```text
Axes
```

[依赖与实现说明](../../faq/bioinformatics.md#optimization)。

## 保存 Matplotlib 图形 {#save-matplotlib-figure}

已有 Matplotlib Figure 可用 `dnakit.visualization.save_figure(figure, target, overwrite=False, tight_layout=True, close=True, **kwargs)` 直接保存。格式由后缀决定，例如 `.svg`、`.png`、`.pdf`；自动创建目录，成功后替换文件，默认整理布局并关闭图形。`dpi` 等额外参数交给 `Figure.savefig`。SVG／PDF 默认省略生成日期；如需自定义 metadata，可显式传入。

```python
import matplotlib.pyplot as plt
from dnakit.visualization import save_figure

figure, axes = plt.subplots()
axes.plot([1, 2, 3], [0.2, 0.4, 0.5])
axes.set(xlabel="Epoch", ylabel="Validation MCC")
save_figure(figure, "output/mcc.svg", overwrite=True)
```

该接口需要 Matplotlib。DNAKit 自身生成的 `SVGArtifact` 继续使用 `save_svg`／`save_image`。
