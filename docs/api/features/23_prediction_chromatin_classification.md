# 表观遗传与开放染色质判别

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) H2A.Z 组蛋白变体区域判别 {#28-h2az}

- **作用：** 判断输入序列是否对应 H2A.Z（H2AFZ）ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h2afz"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
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

## 2) H3K27ac 组蛋白修饰区域判别 {#29-h3k27ac}

- **作用：** 判断输入序列是否对应 H3K27ac ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k27ac"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k27ac",
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

## 3) H3K27me3 组蛋白修饰区域判别 {#30-h3k27me3}

- **作用：** 判断输入序列是否对应 H3K27me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k27me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k27me3",
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

## 4) H3K36me3 组蛋白修饰区域判别 {#31-h3k36me3}

- **作用：** 判断输入序列是否对应 H3K36me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k36me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k36me3",
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

## 5) H3K4me1 组蛋白修饰区域判别 {#32-h3k4me1}

- **作用：** 判断输入序列是否对应 H3K4me1 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k4me1"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me1",
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

## 6) H3K4me2 组蛋白修饰区域判别 {#33-h3k4me2}

- **作用：** 判断输入序列是否对应 H3K4me2 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k4me2"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me2",
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

## 7) H3K4me3 组蛋白修饰区域判别 {#34-h3k4me3}

- **作用：** 判断输入序列是否对应 H3K4me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k4me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me3",
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

## 8) H3K9ac 组蛋白修饰区域判别 {#35-h3k9ac}

- **作用：** 判断输入序列是否对应 H3K9ac ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k9ac"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k9ac",
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

## 9) H3K9me3 组蛋白修饰区域判别 {#36-h3k9me3}

- **作用：** 判断输入序列是否对应 H3K9me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k9me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k9me3",
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
('0', '1')
[0.04569203 0.95430803]
1
```

## 10) H4K20me1 组蛋白修饰区域判别 {#37-h4k20me1}

- **作用：** 判断输入序列是否对应 H4K20me1 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h4k20me1"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h4k20me1",
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

## 11) 人类开放染色质区域判别 {#54}

- **作用：** 判断输入是否属于 Ensembl 人类开放染色质区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_ocr_ensembl"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 80)],
    task="human_ocr_ensembl",
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
