# 可视化

将 DNA 序列、位置、比对、相似度矩阵以及 CSV 中的分类／回归样本分布绘制成可保存的图像。

输入 CSV、指定列名并绘制样本占比图，见下文 [9) 分类与回归样本分布](#viz-009)。

绘图函数返回确定性的正方形 `SVGArtifact`，图内不显示“DNA sequence”等装饰性标题；坐标、高亮标签、图例和碱基数等数据内容仍会保留。配置中的 `title` 仅用于 SVG 无障碍说明，不会绘制在图内。所有图形统一使用 `save_image()` 导出；`image_type` 可选择 `png`、`svg` 或 `jpg`，默认 `png`，目标路径没有扩展名时会自动补充。SVG 不需要额外图形依赖，PNG/JPG 需要安装 `viz` extra；原有 TIFF/PDF 导出继续兼容。下面每段示例都会把图像保存到当前工作目录，“示例结果”展示的是由相同代码生成的真实图像；保存函数默认拒绝覆盖已有文件。

## 1) `VIZ-001` 序列文字图

- **作用：** 把 DNA 碱基、坐标和可选注释绘制成矢量 SVG，返回可嵌入网页或保存的图形，用于快速查看短序列内容。
- **API**：`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`；`config` 使用 `dnakit.visualization.SequencePlotConfig`。
- **输入**：`DNASequence` 或 `DNARecord`；可选每行碱基数、坐标、互补链、行间距、列间距、自定义字符映射和显示上限。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

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

## 2) `VIZ-002` 位置高亮

- **作用**：按坐标、颜色和标签高亮指定碱基或区间，使 motif、突变和其他关注位置在序列图中清晰可见。
- **API**：`dnakit.visualization.Highlight(start[必须], end[必须], label[可选], color[可选], foreground[可选], priority[可选], opacity[可选])`、`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`。
- **输入**：序列和高亮区间；可选标签、颜色和优先级。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

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

## 3) `VIZ-003` 样式控制

- **作用**：通过统一样式配置控制颜色、字体、尺寸、边距、行间距、文字列间距和布局，使不同可视化结果具有一致外观并可复现。
- **API**：`dnakit.visualization.SVGTheme(background[可选], foreground[可选], muted[可选], grid[可选], accent[可选], missing[可选], font_family[可选])`；序列文字图使用 `SequencePlotConfig(font_size[可选], cell_width[可选], line_height[可选], column_spacing[可选], line_spacing[可选], margin[可选], theme[可选])`。
- **输入**：绘图配置；`column_spacing` 控制相邻文字列之间增加的水平像素，`line_spacing` 控制相邻序列行之间增加的垂直像素，二者默认均为 `0`。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

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

## 4) `VIZ-004` Gap 显示

- **作用**：在序列图中用专门符号显示 Gap 的位置、已知或未知长度，同时保持坐标含义，避免把缺失区域误画成真实碱基。
- **API**：`dnakit.Gap(length[必须], kind[可选], crossable[可选], evidence[可选], metadata[可选])`、`dnakit.visualization.plot_sequence(value[必须], highlights[可选], config[可选])`。
- **输入**：包含显式 `Gap` part 的 `DNASequence` 或记录。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

```python
from dnakit import DNASequence, Gap
from dnakit.visualization import plot_sequence, save_image

sequence = DNASequence(["AC", Gap(500), "T", Gap(None), "G"])
artifact = plot_sequence(sequence)
save_image(artifact, "viz-004-gap", image_type="png")
```

- **示例结果（`viz-004-gap.png`）：**

![Gap 显示示例](../../assets/images/visualization/viz-004-gap.png){ width="432" }

## 5) `VIZ-005` 线性序列图

- **作用**：沿线性 DNA 坐标绘制 feature 的起止范围、方向和标签，生成类似基因结构图的 SVG，用于查看多个注释的相对位置。
- **API**：`dnakit.visualization.plot_linear_map(value[必须], width[可选], height[可选], max_features[可选], title[可选], theme[可选])`；画布边长取 `width` 和 `height` 中较大的值，`title` 仅作为 SVG 无障碍说明。
- **输入**：resolved linear `DNASequence` 或带 `DNAFeature` 的 `DNARecord`；输出始终为正方形。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

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

## 6) `VIZ-006` 环状 DNA 图

- **作用**：把环状 DNA 及其 feature 映射到圆周坐标，显示方向、跨原点区间和标签，用于质粒或其他环状分子的结构概览。
- **API**：`dnakit.visualization.plot_circular_map(value[必须], size[可选], max_features[可选], title[可选], theme[可选])`；`size` 是正方形画布边长，`title` 仅作为 SVG 无障碍说明。
- **输入**：`topology="circular"` 的 resolved 序列或记录。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

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

## 7) `VIZ-007` 自定义字符表示

- **作用**：把 A、T、C、G 或 IUPAC 符号替换为调用方指定的显示字符，例如 A 显示为 `*`、T 显示为 `-`、C 显示为 `+`、G 显示为 `]`；坐标、高亮、互补关系和碱基颜色仍按原始序列计算。
- **API**：`dnakit.visualization.SequencePlotConfig(symbol_map[可选])`；通过 `plot_sequence(value, config=config)` 使用。
- **输入**：`symbol_map` 的键是大写 DNA/IUPAC 符号，值是一个可见字符；未配置的符号保持原样，默认映射为空。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

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

## 8) `VIZ-008` Alignment 图

- **作用**：将成对比对的两条序列按列绘制，区分匹配、错配、插入和删除，并显示坐标，便于人工检查差异位置。
- **API**：`dnakit.visualization.plot_alignment(result[必须], columns_per_line[可选], max_columns[可选], theme[可选])`。
- **输入**：预先计算的 `AlignmentResult`；可选每行列数和上限，输出始终为正方形。
- **导出**：`save_image(artifact, target, image_type[可选])`；可选 `png`、`svg`、`jpg`，默认 `png`。
- **示例代码**：

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

## 9) `VIZ-009` 分类与回归样本分布 {#viz-009}

读取 CSV 中的一列或多列，为每列单独绘制样本占比图。**柱高 = 该类别或区间内的样本数 ÷ 该列有效样本总数**，纵轴以百分比显示，每张图的柱高之和为 100%（显示值四舍五入后可能略有偏差）。这是数据分布，不是模型对单个样本的预测置信度。

| 模式 | 输入列 | 柱子含义 |
| --- | --- | --- |
| `classification` | `0`、`1`、`2` 等类别标签，也支持字符串 | 每个类别一根柱；类别数量自动扩展 |
| `regression` | 连续数值 | 将数值分箱，每个区间一根柱；默认 30 个等宽区间 |

### CSV 示例与统计规则

完整的 [示例 CSV](../../assets/data/distribution.csv) 位于仓库的 `docs/assets/data/distribution.csv`，包含 200 条人工生成的演示数据：`label` 列的 0、1、2 类分别有 40、100、60 条；`score`、`score_alt` 是连续数值。CSV 首行必须是列名，绘图时按列名精确选择。

```text
sample_id,label,score,score_alt
sample_001,0,2.6656,6.5507
sample_002,0,6.0421,6.6343
```

- 每列独立统计，互不混合；多列在一次调用中使用相同模式和样式，分类与回归混合时分别调用。
- 空单元格、纯空白字符串和 Python 输入中的 `None` 跳过；缺失数量记录在结果中，CLI 也会报告。整列为空时报错。
- 回归的非数字、`NaN`、`inf` 报错，不自动当作缺失值；分类保留标签文本，`0` 与 `0.0` 是不同类别，首尾空白会去除。
- 分类默认按类别首次出现顺序排列，可指定顺序并增加未出现的类别，后者占比为 0%；指定类别必须包含全部已观察到的类别。
- 回归区间为左闭右开 `[a, b)`，最后一个区间包含右端点。常数列自动画一个占比 100% 的区间。
- 自定义分箱边界必须严格递增并覆盖全部有效数值，超出范围时报错；即使区间不等宽，柱高仍是样本占比，**不是概率密度**。
- 默认最多读取 1,000,000 行、绘制 1,000 个类别／区间，超限时报错，不截断统计。刻度过多时会减少文字标签，所有柱子仍保留。

### 命令行：分类占比

在项目根目录、已安装 DNAKit 的环境中执行。若尚未安装，先运行 `python -m pip install -e ".[viz]"`；SVG 导出仅需基础安装，PNG/JPG/PDF 等使用 `viz` 依赖。

作用：统计 `label` 列的类别占比，按 0、1、2 排列，柱上显示百分比，并保存为 PNG。

```bash
python -m dnakit.cli plot-distribution \
  --input docs/assets/data/distribution.csv \
  --column label \
  --task classification \
  --output outputs/label_distribution.png \
  --category 0 \
  --category 1 \
  --category 2 \
  --show-values \
  --font-family sans-serif \
  --font-size 18 \
  --color '#2563eb' \
  --dpi 150
```

示例结果：三根柱分别为 20%、50%、30%。

![分类样本占比](../../assets/images/visualization/viz-009-classification.png){ width="480" }

### 命令行：回归数值分布

作用：对 `score` 列按 30 个等宽区间统计占比，使用橙色小柱保存为 PNG。

```bash
python -m dnakit.cli plot-distribution \
  --input docs/assets/data/distribution.csv \
  --column score \
  --task regression \
  --output outputs/score_distribution.png \
  --bins 30 \
  --bar-width 0.95 \
  --font-size 18 \
  --color '#ea580c' \
  --dpi 150
```

![回归样本占比直方图](../../assets/images/visualization/viz-009-regression.png){ width="480" }

作用：一次读取两列，各自统计并保存为指定的 SVG 文件。`--column` 与 `--output` 按出现顺序一一对应。

```bash
python -m dnakit.cli plot-distribution \
  --input docs/assets/data/distribution.csv \
  --column score \
  --column score_alt \
  --task regression \
  --output outputs/score.svg \
  --output outputs/score_alt.svg \
  --bins 25
```

自定义区间时使用 `--bin-edges '0,2,4,6,8,10'`，优先于 `--bins`。边界从负数开始时使用 `--bin-edges=-10,-5,0,5,10,15`，避免被命令行解析为选项。示例两列分别自动确定范围；要比较相同区间，应提供共同的 `--bin-edges`。

### Python API 与样式

`plot_csv_distributions(source, columns=..., task=..., config=...)` 返回 `{列名: SVGArtifact}`；已经有列表时，可调用 `plot_distribution(values, task=..., config=...)`。使用既有 `save_image()` 指定文件名和格式。

```python
from dnakit.visualization import (
    DistributionPlotConfig,
    ImageExportConfig,
    SVGTheme,
    plot_csv_distributions,
    save_image,
)

plots = plot_csv_distributions(
    "docs/assets/data/distribution.csv",
    columns=["label"],
    task="classification",
    config=DistributionPlotConfig(
        categories=("0", "1", "2", "3"),  # 未出现的第 3 类占比为 0%
        size=800,
        font_size=20,
        tick_font_size=16,
        value_font_size=14,
        show_values=True,
        bar_width=0.8,
        x_label="Class",
        theme=SVGTheme(font_family="sans-serif", accent="#2563eb"),
    ),
)
artifact = plots["label"]
print(artifact.metadata["counts"])       # (40, 100, 60, 0)
print(artifact.metadata["proportions"])  # (0.2, 0.5, 0.3, 0.0)
save_image(
    artifact,
    "outputs/custom_distribution.png",
    config=ImageExportConfig(dpi=150, create_parents=True),
)
```

| 配置／CLI 参数 | 默认值 | 作用 |
| --- | --- | --- |
| `bins` / `--bins` | `30` | 回归等宽分箱数量 |
| `bin_edges` / `--bin-edges` | `None` | 回归自定义边界，覆盖 `bins` |
| `categories` / 重复 `--category` | `None` | 类别显示顺序及可选的零计数类别 |
| `size` / `--size` | `720` | 正方形 SVG 画布边长，至少 320 |
| `font_size` / `--font-size` | `16` | 坐标轴名称字号 |
| `tick_font_size` / `--tick-font-size` | `12` | 坐标刻度字号 |
| `value_font_size` / `--value-font-size` | `11` | 柱顶占比字号 |
| `theme.font_family` / `--font-family` | API `monospace`；CLI `sans-serif` | 字体名称，需系统已安装；中文可使用 `Noto Sans CJK SC` |
| `theme.accent` / `--color` | `#2563eb` | 柱子颜色 |
| `theme.background` / `--background` | `#ffffff` | 背景颜色；API 的 `SVGTheme` 还可设置文字、网格颜色 |
| `bar_width` / `--bar-width` | `0.9` | 柱宽占类别槽或分箱宽度的比例，`1` 表示无间隙 |
| `label_rotation` / `--label-rotation` | `0` | 横轴刻度旋转角度，可选 `0`、`45`、`90` |
| `max_tick_labels` / `--max-tick-labels` | `12` | 控制横轴标签密度；回归额外保留最右边界 |
| `show_values` / `--show-values` | `False` | 显示柱顶占比；密集直方图建议保持关闭 |
| `show_grid` / `--grid`、`--no-grid` | `True` | 显示水平网格 |
| `x_label` / `--x-label` | CSV 列名 | 横轴名称 |
| `y_label` / `--y-label` | `Sample proportion (%)` | 纵轴名称；数据始终按占比计算 |
| `ImageExportConfig.dpi` / `--dpi` | `600` | 位图导出分辨率；默认像素边长约为 `size × dpi / 96` |
| `encoding` / `--encoding` | `utf-8-sig` | CSV 编码，默认兼容 UTF-8 BOM |
| `delimiter` / `--delimiter` | `,` | CSV 单字符分隔符 |
| `max_values` / `--max-values` | `1000000` | 输入样本行数上限，包含空值行 |

输出支持 PNG、SVG、JPG、TIFF、PDF；没有扩展名时默认补 `.png`。CLI 自动创建输出父目录，默认拒绝覆盖，重复执行同一路径时可加 `--overwrite`。与现有可视化保持一致，`title` 仅是 SVG 无障碍说明，不在图内绘制标题。
