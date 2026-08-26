# 可视化

将 DNA 序列、位置、比对和相似度矩阵绘制成可保存的图像。

绘图函数返回确定性的 `SVGArtifact`；SVG 不需要额外图形依赖，PNG/TIFF/PDF 需要安装 `viz` extra。下面每段示例都会把图像保存到当前工作目录，“示例结果”展示的是由相同代码生成的真实图像；保存函数默认拒绝覆盖已有文件。

## 1) `VIZ-001` 序列文字图

- **作用：** 把 DNA 碱基、坐标和可选注释绘制成矢量 SVG，返回可嵌入网页或保存的图形，用于快速查看短序列内容。
- **API**：`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`；`config` 使用 `dnakit.visualization.SequencePlotConfig`。
- **输入**：`DNASequence` 或 `DNARecord`；可选每行碱基数、坐标、互补链和显示上限。
- **示例代码**：

```python
from dnakit import DNASequence
from dnakit.visualization import SequencePlotConfig, plot_sequence, save_svg

artifact = plot_sequence(
    DNASequence("ACGTACGT"),
    config=SequencePlotConfig(bases_per_line=4),
)
save_svg(artifact, "viz-001-sequence.svg")
```

- **示例结果（`viz-001-sequence.svg`）：**

![序列文字图示例](../../assets/images/visualization/viz-001-sequence.svg)

- **限制**：超过上限默认报错；只有显式选择截断策略才会截断显示。

## 2) `VIZ-002` 位置高亮

- **作用**：按坐标、颜色和标签高亮指定碱基或区间，使 motif、突变和其他关注位置在序列图中清晰可见。
- **API**：`dnakit.visualization.Highlight(start[必须], end[必须], label[可选], color[可选], foreground[可选], priority[可选], opacity[可选])`、`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`。
- **输入**：序列和高亮区间；可选标签、颜色和优先级。
- **示例代码**：

```python
from dnakit import DNASequence
from dnakit.visualization import Highlight, plot_sequence, save_svg

artifact = plot_sequence(
    DNASequence("ACGTACGT"),
    highlights=[Highlight(1, 4, label="target", color="#ffcc00")],
)
save_svg(artifact, "viz-002-highlight.svg")
```

- **示例结果（`viz-002-highlight.svg`）：**

![位置高亮示例](../../assets/images/visualization/viz-002-highlight.svg)

- **限制**：坐标相对于 `DNASequence.symbols`，不把显式 Gap 长度计入索引。

## 3) `VIZ-003` 样式控制

- **作用**：通过统一样式配置控制颜色、字体、尺寸、边距和布局，使不同可视化结果具有一致外观并可复现。
- **API**：`dnakit.visualization.SVGTheme(background[可选], foreground[可选], muted[可选], grid[可选], accent[可选], missing[可选], font_family[可选])`；各绘图函数通过其可选 `config` 或 `theme` 参数接收样式。
- **输入**：绘图配置；可选颜色、字体和尺寸。
- **示例代码**：

```python
from dnakit import DNASequence
from dnakit.visualization import SVGTheme, SequencePlotConfig, plot_sequence, save_svg

theme = SVGTheme(background="#111827", foreground="#f9fafb", accent="#22c55e")
artifact = plot_sequence(
    DNASequence("ACGT"),
    config=SequencePlotConfig(font_size=18, theme=theme),
)
save_svg(artifact, "viz-003-theme.svg")
```

- **示例结果（`viz-003-theme.svg`）：**

![深色样式控制示例](../../assets/images/visualization/viz-003-theme.svg)

- **限制**：颜色和字体名称接受严格格式校验。

## 4) `VIZ-004` Gap 显示

- **作用**：在序列图中用专门符号显示 Gap 的位置、已知或未知长度，同时保持坐标含义，避免把缺失区域误画成真实碱基。
- **API**：`dnakit.Gap(length[必须], kind[可选], crossable[可选], evidence[可选], metadata[可选])`、`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`。
- **输入**：包含显式 `Gap` part 的 `DNASequence` 或记录。
- **示例代码**：

```python
from dnakit import DNASequence, Gap
from dnakit.visualization import plot_sequence, save_svg

sequence = DNASequence(["AC", Gap(500), "T", Gap(None), "G"])
artifact = plot_sequence(sequence)
save_svg(artifact, "viz-004-gap.svg")
```

- **示例结果（`viz-004-gap.svg`）：**

![Gap 显示示例](../../assets/images/visualization/viz-004-gap.svg)

- **限制**：Gap 不会为绘图而展开成大量字符；未知 Gap 后坐标显示为未知。

## 5) `VIZ-005` 线性序列图

- **作用**：沿线性 DNA 坐标绘制 feature 的起止范围、方向和标签，生成类似基因结构图的 SVG，用于查看多个注释的相对位置。
- **API**：`dnakit.visualization.plot_linear_map(value[必须], width[可选], height[可选], max_features[可选], title[可选], theme[可选])`。
- **输入**：resolved linear `DNASequence` 或带 `DNAFeature` 的 `DNARecord`。
- **示例代码**：

```python
from dnakit import DNAFeature, DNARecord, DNASequence, Interval
from dnakit.visualization import plot_linear_map, save_svg

record = DNARecord(
    DNASequence("ACGTACGT"),
    "linear",
    features=[DNAFeature("gene", Interval(1, 6))],
)
artifact = plot_linear_map(record)
save_svg(artifact, "viz-005-linear-map.svg")
```

- **示例结果（`viz-005-linear-map.svg`）：**

![线性序列图示例](../../assets/images/visualization/viz-005-linear-map.svg)

- **限制**：要求可解析的线性 coordinate span；未解析 feature 不会被猜测绘制。

## 6) `VIZ-006` 环状 DNA 图

- **作用**：把环状 DNA 及其 feature 映射到圆周坐标，显示方向、跨原点区间和标签，用于质粒或其他环状分子的结构概览。
- **API**：`dnakit.visualization.plot_circular_map(value[必须], size[可选], max_features[可选], title[可选], theme[可选])`。
- **输入**：`topology="circular"` 的 resolved 序列或记录。
- **示例代码**：

```python
from dnakit import DNAFeature, DNARecord, DNASequence, Interval
from dnakit.visualization import plot_circular_map, save_svg

record = DNARecord(
    DNASequence("ACGTACGT", topology="circular"),
    "plasmid",
    features=[DNAFeature("gene", Interval(1, 6))],
)
artifact = plot_circular_map(record)
save_svg(artifact, "viz-006-circular-map.svg")
```

- **示例结果（`viz-006-circular-map.svg`）：**

![环状 DNA 图示例](../../assets/images/visualization/viz-006-circular-map.svg)

- **限制**：不会自动选择新起点或计算未解析 feature。

## 7) `VIZ-008` Alignment 图

- **作用**：将成对比对的两条序列按列绘制，区分匹配、错配、插入和删除，并显示坐标，便于人工检查差异位置。
- **API**：`dnakit.visualization.plot_alignment(result[必须], columns_per_line[可选], max_columns[可选], theme[可选])`。
- **输入**：预先计算的 `AlignmentResult`；可选每行列数和上限。
- **示例代码**：

```python
from dnakit import DNASequence
from dnakit.alignment import align_pairwise
from dnakit.visualization import plot_alignment, save_svg

alignment = align_pairwise(DNASequence("ACGT"), DNASequence("AGT"))
artifact = plot_alignment(alignment, columns_per_line=4)
save_svg(artifact, "viz-008-alignment.svg")
```

- **示例结果（`viz-008-alignment.svg`）：**

![Alignment 图示例](../../assets/images/visualization/viz-008-alignment.svg)

- **限制**：函数只消费既有结果，不负责选择比对参数。
