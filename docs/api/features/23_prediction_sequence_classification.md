# 物种与分子功能分类

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) SupKTax 分类 {#19-supktax}

- **作用：** 使用官方 SupKTax checkpoint 对 gene 序列进行分类。
- **预测方法：** 调用 LucaOneTasks 官方 `SupKTax` 多分类 checkpoint 直接推理，并返回概率最高的 `top_k` 个 `label.txt` 标签。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="supktax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`。结果标签以该 checkpoint 的 `label.txt` 为准。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="supktax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个官方标签>,)
```

## 2) GenusTax 分类 {#20-genustax}

- **作用：** 使用官方 GenusTax checkpoint 对 gene 序列进行分类。
- **预测方法：** 调用 LucaOneTasks 官方 `GenusTax` 多分类 checkpoint 直接推理，并按概率返回 Top-k 标签。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="genustax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`。LucaOneTasks README 对 SupKTax/GenusTax 的自然语言层级存在对调，因此解释结果时以实际 `label.txt` 为准。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="genustax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["predicted_label"])
```

- **示例结果（结构示意）：**

```text
(5,)
<概率最高的官方标签>
```

## 3) 物种分类 {#21}

- **作用：** 预测 gene 序列对应的物种类别。
- **预测方法：** 调用 LucaOneTasks 官方 `SpeciesTax` 多分类 checkpoint 直接推理，返回 Top-k 物种标签和对应概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="speciestax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`；可用 `top_k` 设置返回标签数。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="speciestax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个物种标签>,)
```

## 4) 原核蛋白亚细胞定位 {#22}

- **作用：** 预测原核蛋白的亚细胞定位类别。
- **预测方法：** 调用 LucaOneTasks 官方 `ProtLoc` 多分类 checkpoint 直接推理，返回 Top-k 定位标签及概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_location", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="protein"` 的 `BiologicalSequence`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [
        BiologicalSequence(
            "protein-1",
            "MAMAPRTEINSTRING" * 10,
            "protein",
        )
    ],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_location",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个亚细胞定位标签>,)
```

## 5) ncRNA 家族 {#24-ncrna}

- **作用：** 预测非编码 RNA 序列所属的家族。
- **预测方法：** 调用 LucaOneTasks 官方 `ncRNAFam` 多分类 checkpoint 直接推理，返回 Top-k 家族标签和概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="ncrna_family", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`；DNA 或 RNA 字母均可，`U` 由官方任务流程处理。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("ncrna-1", "AUGCUA" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="ncrna_family",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个 ncRNA 家族标签>,)
```

## 6) 人类与线虫序列判别 {#47}

- **作用：** 判断输入基因组片段来自人类还是秀丽隐杆线虫。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="demo_human_or_worm"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 50)],
    task="demo_human_or_worm",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```
