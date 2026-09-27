# 启动子与增强子判别

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) 增强子判别 {#38}

- **作用：** 判断输入序列属于增强子还是非增强子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="enhancers"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="enhancers",
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

## 2) 增强子类型分类 {#39}

- **作用：** 在非增强子、组织特异增强子和组织不变增强子之间进行三分类。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 3 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="enhancers_types"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="enhancers_types",
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
(<类别0概率>, <类别1概率>, <类别2概率>)
<概率最高的 checkpoint 标签>
```

## 3) 启动子判别 {#40}

- **作用：** 判断输入序列属于启动子还是非启动子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="promoter_all"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_all",
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

## 4) 非 TATA 启动子判别 {#41-tata}

- **作用：** 判断输入序列属于非 TATA 启动子还是负样本。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="promoter_no_tata"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_no_tata",
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

## 5) TATA 启动子判别 {#42-tata}

- **作用：** 判断输入序列属于 TATA 启动子还是负样本。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="promoter_tata"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_tata",
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

## 6) 果蝇增强子判别 {#48}

- **作用：** 判断输入是否属于 Stark 果蝇增强子集合。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="drosophila_enhancers_stark"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 500)],
    task="drosophila_enhancers_stark",
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

## 7) 小鼠增强子判别 {#49}

- **作用：** 判断输入是否属于 Ensembl 小鼠增强子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="dummy_mouse_enhancers_ensembl"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 500)],
    task="dummy_mouse_enhancers_ensembl",
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

## 8) 人类 Cohn 增强子判别 {#50-cohn}

- **作用：** 判断输入是否属于 Cohn 人类增强子集合。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_enhancers_cohn"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 125)],
    task="human_enhancers_cohn",
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

## 9) 人类 Ensembl 增强子判别 {#51-ensembl}

- **作用：** 判断输入是否属于 Ensembl 人类增强子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_enhancers_ensembl"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="human_enhancers_ensembl",
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

## 10) 人类调控元件类型分类 {#52}

- **作用：** 在人类增强子、启动子和开放染色质区域之间进行三分类。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 3 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_ensembl_regulatory"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="human_ensembl_regulatory",
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
(<类别0概率>, <类别1概率>, <类别2概率>)
<概率最高的 checkpoint 标签>
```

## 11) 人类非 TATA 启动子判别 {#53-tata}

- **作用：** 判断输入是否属于人类非 TATA 启动子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_nontata_promoters"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", ("ACGT" * 63)[:251])],
    task="human_nontata_promoters",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（已有本地 checkpoint 的 CUDA smoke 记录）：** 依次输出类别名称、概率和预测标签；数值随序列与 checkpoint 而变化。

```text
('negative', 'positive')
[3.8891116e-08 1.0000000e+00]
positive
```
