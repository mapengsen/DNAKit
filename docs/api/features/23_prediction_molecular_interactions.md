# 分子关系与互作预测

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) 核酸–蛋白中心法则关系 {#18}

- **作用：** 预测一条 gene 序列与一条 protein 序列是否具有 CentralDogma 关系。
- **预测方法：** 调用 LucaOneTasks 官方 `CentralDogma` 二分类 checkpoint 和 `src/predict_v1.py` 直接推理，返回正类概率及官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="central_dogma", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；第一条必须为 `sequence_type="gene"`，第二条必须为 `sequence_type="protein"`。`model_source_path` 必须指向 LucaOneTasks 官方源码目录。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "central-dogma-1",
    BiologicalSequence("gene-1", "ATGGCC" * 100, "gene"),
    BiologicalSequence("protein-1", "MAMAPRTEINSTRING" * 10, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="central_dogma",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```

## 2) 甲型流感抗原关系 {#25}

- **作用：** 预测两条甲型流感 gene 序列之间的抗原关系类别。
- **预测方法：** 调用 LucaOneTasks 官方 `InfA` 二分类 checkpoint 直接推理，返回正类概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="influenza_antigenicity", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；两条序列都必须是 `sequence_type="gene"`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "influenza-pair-1",
    BiologicalSequence("strain-a", "ATGGCC" * 100, "gene"),
    BiologicalSequence("strain-b", "ATGGCT" * 100, "gene"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="influenza_antigenicity",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```

## 3) 蛋白–蛋白互作 {#26}

- **作用：** 预测两条蛋白序列是否具有蛋白–蛋白相互作用。
- **预测方法：** 调用 LucaOneTasks 官方 `PPI` 二分类 checkpoint 直接推理，返回互作概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_interaction", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；两条序列都必须是 `sequence_type="protein"`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "protein-pair-1",
    BiologicalSequence("protein-a", "MAMAPRTEINSTRING" * 10, "protein"),
    BiologicalSequence("protein-b", "MKWVTFISLLFLFSSAYSR" * 8, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_interaction",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```

## 4) ncRNA–蛋白互作 {#27-ncrna}

- **作用：** 预测一条 ncRNA 与一条蛋白序列是否具有相互作用。
- **预测方法：** 调用 LucaOneTasks 官方 `ncRPI` 二分类 checkpoint 直接推理，返回互作概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="ncrna_protein_interaction", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；第一条必须是 `sequence_type="gene"`，第二条必须是 `sequence_type="protein"`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "ncrna-protein-1",
    BiologicalSequence("ncrna-1", "AUGCUA" * 100, "gene"),
    BiologicalSequence("protein-1", "MKWVTFISLLFLFSSAYSR" * 8, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="ncrna_protein_interaction",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```
