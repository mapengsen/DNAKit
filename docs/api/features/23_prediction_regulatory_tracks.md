# 表达、调控与接触图预测

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) RNA-seq 覆盖与表达 {#1-rna-seq}

- **作用：** 根据 DNA 序列预测 RNA-seq 覆盖和基因表达轨道。
- **预测方法：** 使用 AlphaGenome all-folds 官方权重直接推理 `RNA_SEQ` 输出；短序列会居中补 `N` 到官方支持的上下文长度，可用 `ontology_terms` 限定组织或细胞类型。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="rna_seq")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也可改用 `VariantContext` 和 `predict_variant_effects()`，得到 REF、ALT、ALT−REF 三组轨道。支持人或小鼠。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="rna_seq",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <RNA-seq 轨道数>)
('position', 'track')
```

## 2) CAGE 转录起始与表达 {#2-cage}

- **作用：** 预测 CAGE 信号，用于描述加帽转录本的转录起始位置和表达强度。
- **预测方法：** 使用 AlphaGenome 官方 `CAGE` 输出头直接生成位置 × 轨道矩阵；`ontology_terms` 可筛选目标组织或细胞类型。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="cage")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext` 形式的单 SNV 上下文。序列最长 1,048,576 bp，支持人或小鼠。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="cage",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <CAGE 轨道数>)
('position', 'track')
```

## 3) PRO-cap 新生转录起始 {#3-pro-cap}

- **作用：** 预测 PRO-cap 新生 RNA 5′ 端信号，用于定位活跃转录起始事件。
- **预测方法：** 使用 AlphaGenome 官方 `PROCAP` 输出头直接推理；输出分辨率记录在 `metadata["resolution_bp"]` 中。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="procap")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可指定 `organism` 和 `ontology_terms`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="procap",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["resolution_bp"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <PRO-cap 轨道数>)
<每个输出位置对应的 bp 数>
```

## 4) ATAC-seq 染色质可及性 {#4-atac-seq}

- **作用：** 预测 ATAC-seq 开放染色质轨道，用于识别可能处于可及状态的调控区域。
- **预测方法：** 使用 AlphaGenome 官方 `ATAC` 输出头直接推理，并按位置返回不同组织或细胞类型的可及性信号。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="atac")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持只包含一个碱基变化的 `VariantContext`。可通过 `ontology_terms` 限制轨道。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="atac",
        ontology_terms=("UBERON:0002048",),
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names[:2])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <筛选后的 ATAC-seq 轨道数>)
(<前两个官方轨道名称>,)
```

## 5) DNase-seq 染色质可及性 {#5-dnase-seq}

- **作用：** 预测 DNase-seq 超敏信号，用于描述染色质可及性。
- **预测方法：** 使用 AlphaGenome 官方 `DNASE` 输出头直接推理；短输入会居中补齐，补齐长度保存在结果元数据中。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="dnase")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可选 `organism="human"` 或 `"mouse"`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="dnase",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["input_length"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <DNase-seq 轨道数>)
16384
```

## 6) 组蛋白修饰 {#6}

- **作用：** 预测组蛋白 ChIP-seq 轨道，用于表征不同组蛋白修饰在序列区域中的信号。
- **预测方法：** 使用 AlphaGenome 官方 `CHIP_HISTONE` 输出头直接推理，返回位置 × 组蛋白轨道矩阵。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="chip_histone")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可按物种和 ontology term 选择输出。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="chip_histone",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <组蛋白 ChIP-seq 轨道数>)
('position', 'track')
```

## 7) 转录因子结合 {#7}

- **作用：** 预测转录因子 ChIP-seq 结合轨道，用于定位可能的转录因子结合信号。
- **预测方法：** 使用 AlphaGenome 官方 `CHIP_TF` 输出头直接推理，返回位置 × 转录因子轨道矩阵。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="chip_tf")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可用 `ontology_terms` 限定组织或细胞类型。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="chip_tf",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names[:2])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <TF ChIP-seq 轨道数>)
(<前两个官方轨道名称>,)
```

## 8) 三维 DNA 接触图 {#11-dna}

- **作用：** 预测染色质区域之间的三维 DNA 接触概率图。
- **预测方法：** 使用 AlphaGenome 官方 `CONTACT_MAPS` 输出头直接推理，输出两个位置轴和一个轨道轴。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="contact_maps")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`，变异结果在最前面增加 REF、ALT、ALT−REF 轴。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="contact_maps",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<bin 数>, <bin 数>, <接触图轨道数>)
('position_1', 'position_2', 'track')
```

## 9) 人类调控轨道 {#12}

- **作用：** 预测人类基因组的 5,313 条 CAGE、DNase/ATAC、转录因子 ChIP 和组蛋白 ChIP 轨道。
- **预测方法：** 使用公开的 Enformer 预训练权重直接推理；默认把输入居中补 `N` 到 196,608 bp，并输出 896 个 128 bp bin。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="enformer", task="human_tracks")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也可用 `VariantContext` 和 `predict_variant_effects()`，获得 REF、ALT、ALT−REF。`max_length` 必须是 114,688–393,216 范围内的 128 倍数。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="enformer",
        task="human_tracks",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["bin_size_bp"])
```

- **示例结果（默认上下文）：**

```text
(896, 5313)
128
```

## 10) 小鼠调控轨道 {#13}

- **作用：** 预测小鼠基因组的 1,643 条 CAGE、DNase/ATAC、转录因子 ChIP 和组蛋白 ChIP 轨道。
- **预测方法：** 使用 Enformer 的小鼠输出头直接推理；默认上下文为 196,608 bp，输出为 128 bp 分箱轨道。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="enformer", task="mouse_tracks")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持只含一个 SNV 的 `VariantContext`。输入不足上下文长度时自动在两端居中补 `N`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("mouse-region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="enformer",
        task="mouse_tracks",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["organism"])
```

- **示例结果（默认上下文）：**

```text
(896, 1643)
mouse
```

<!-- BEGIN GENERATED API: model-tracks -->

<span id="model-tracks"></span>

## 11) 预测基因组区间的多模态调控轨道 {#dnakit.databases.alpha_genome_predict_interval}

- **作用：** 预测基因组区间的多模态调控轨道。
- **处理方法：** 通过 `AlphaGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.alpha_genome_predict_interval(*, chromosome, start, end, output_types=None, ontology_terms=None, organism=None, sequence_length=None, verbose=True)`
- **输入：** 必须提供 `chromosome`, `start`, `end`。 需配置：`ALPHA_GENOME_API_KEY`.
- **输出：** 返回字典，包含 `interval`, `tracks` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.alpha_genome_predict_interval(
    chromosome='chr19',
    start=40991281,
    end=41018398,
    output_types=['RNA_SEQ', 'ATAC'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "interval": "<string>",
  "tracks": [
    "<object>"
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_regulatory_tracks.md#dnakit.databases.alpha_genome_predict_interval-parameters)

<!-- END GENERATED API: model-tracks -->
