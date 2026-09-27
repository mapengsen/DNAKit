# 剪接与基因组区域预测

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) 剪接供体/受体位点 {#8}

- **作用：** 预测序列各位置的剪接供体和剪接受体相关轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_SITES` 输出头直接推理，并保留输出分辨率及轨道名称。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_sites")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`，变异结果额外包含 ALT−REF 差值轴。
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
        task="splice_sites",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置数>, <供体/受体轨道数>)
('position', 'track')
```

## 2) 剪接位点使用比例 {#9}

- **作用：** 预测候选剪接位点的使用信号或使用比例轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_SITE_USAGE` 输出头直接推理；具体轨道身份由 `output_names` 给出。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_site_usage")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可指定物种、ontology term 和官方上下文长度。
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
        task="splice_site_usage",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["resolution_bp"])
```

- **示例结果（结构示意）：**

```text
(<位置数>, <剪接位点使用轨道数>)
<每个输出位置对应的 bp 数>
```

## 3) 剪接连接计数 {#10}

- **作用：** 预测候选 splice junction 的连接计数轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_JUNCTIONS` 输出头直接推理；结果元数据中的 `junctions` 保存候选连接标识。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_junctions")[必须], backend[可选])`
- **输入：** 只接受 `BiologicalSequence`。该任务的候选连接集合可能随等位基因变化，因此当前不接受 `VariantContext`。
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
        task="splice_junctions",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
print(len(output.metadata.get("junctions", ())))
```

- **示例结果（结构示意）：**

```text
(<候选剪接连接数>, <轨道数>)
('junction', 'track')
<候选剪接连接数>
```

## 4) 14 类单碱基基因组分割 {#14-14}

- **作用：** 为人 DNA 序列中的每个碱基预测 14 类基因和调控元件概率。
- **预测方法：** 使用带训练完成分割头的 `InstaDeepAI/segment_nt`，即 NT-v2 骨干网络加 U-Net 分割头；对二分类 logits 做 softmax 后返回“该特征存在”的概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="segmentnt", task="genomic_segmentation", allow_remote_code=True)[必须], backend[可选])`
- **输入：** 最长 30,000 bp 的人 `BiologicalSequence`。14 个输出类别见下表。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 1000)],
    config=PropertyPredictionConfig(
        model="segmentnt",
        task="genomic_segmentation",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(len(output.output_names))
```

- **示例结果：**

```text
(4000, 14)
14
```

以下 14 个类别共用上述预测接口，每类对应 `output.values` 的一列，实际列名与顺序以 `output.output_names` 为准。

| 类别 | 含义 |
| --- | --- |
| <span id="141-protein-coding-gene"></span>蛋白质编码基因（protein-coding gene） | 预测每个碱基属于蛋白质编码基因区域的概率。 |
| <span id="142-rnalncrna"></span>长链非编码 RNA（lncRNA） | 预测每个碱基属于长链非编码 RNA 基因区域的概率。 |
| <span id="143-exon"></span>外显子（exon） | 预测每个碱基属于外显子区域的概率。 |
| <span id="144-intron"></span>内含子（intron） | 预测每个碱基属于内含子区域的概率。 |
| <span id="145-splice-donor"></span>剪接供体位点（splice donor） | 预测每个碱基属于剪接供体位点的概率。 |
| <span id="146-splice-acceptor"></span>剪接受体位点（splice acceptor） | 预测每个碱基属于剪接受体位点的概率。 |
| <span id="147-5-5-utr"></span>5′ 非翻译区（5′ UTR） | 预测每个碱基属于 5′ 非翻译区的概率。 |
| <span id="148-3-3-utr"></span>3′ 非翻译区（3′ UTR） | 预测每个碱基属于 3′ 非翻译区的概率。 |
| <span id="149-ctcf-ctcf-bound"></span>CTCF 结合区域（CTCF-bound） | 预测每个碱基属于 CTCF 结合区域的概率。 |
| <span id="1410-polya-signal"></span>多聚腺苷酸化信号（polyA signal） | 预测每个碱基属于多聚腺苷酸化信号的概率。 |
| <span id="1411-tissue-specific-enhancer"></span>组织特异性增强子（tissue-specific enhancer） | 预测每个碱基属于组织特异性增强子区域的概率。 |
| <span id="1412-tissue-invariant-enhancer"></span>组织不变性增强子（tissue-invariant enhancer） | 预测每个碱基属于组织不变性增强子区域的概率。 |
| <span id="1413-tissue-specific-promoter"></span>组织特异性启动子（tissue-specific promoter） | 预测每个碱基属于组织特异性启动子区域的概率。 |
| <span id="1414-tissue-invariant-promoter"></span>组织不变性启动子（tissue-invariant promoter） | 预测每个碱基属于组织不变性启动子区域的概率。 |

## 5) 外显子概率 {#16}

- **作用：** 根据一个候选位点两侧的正向和反向上下文预测该位点属于 exon 的概率。
- **预测方法：** 从 Evo 2 `blocks.26` 读取正向、反向上下文的末 token 4,096 维 embedding，拼接为 8,192 维，再通过官方已经训练好的 MLP 头输出概率；不会重新训练分类器。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="evo2", task="exon_probability", allow_remote_code=True)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`，顺序必须是 forward context、reverse context；两条都使用 `sequence_type="gene"`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

contexts = BiologicalSequencePair(
    "site-1",
    BiologicalSequence("forward-1", "ACGT" * 256),
    BiologicalSequence("reverse-1", "ACGT" * 256),
)
result = predict_pair_properties(
    [contexts],
    config=PropertyPredictionConfig(
        model="evo2",
        task="exon_probability",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签取决于概率和阈值）：**

```text
('exon_probability',)
exon 或 non_exon
```

## 6) 剪接受体位点判别 {#43}

- **作用：** 判断序列中央上下文是否包含剪接受体位点。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="splice_sites_acceptors"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_acceptors",
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

## 7) 剪接位点类型分类 {#44}

- **作用：** 在无剪接位点、剪接受体和剪接供体之间进行三分类。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 3 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="splice_sites_all"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_all",
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

## 8) 剪接供体位点判别 {#45}

- **作用：** 判断序列中央上下文是否包含剪接供体位点。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="splice_sites_donors"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_donors",
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

## 9) 编码区与基因间区判别 {#46}

- **作用：** 判断输入来自蛋白编码区还是基因间区。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="demo_coding_vs_intergenomic_seqs"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 50)],
    task="demo_coding_vs_intergenomic_seqs",
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
