# 深度学习性质预测

本页只收录**官方已经发布预训练任务头、可以直接推理**的功能。DNAKit 只加载现成权重，不会在本地重新训练或微调模型。原有的[神经网络表征](08_fingerprints.md)仍只负责提取 embedding。

下面每个功能都单独说明作用、预测方法、API、输入与输出。示例结果中的占位符只表示输出结构；具体数值、标签数和轨道数取决于输入、物种、组织筛选条件及实际 checkpoint。

本章涉及的模型论文、官方仓库和权重来源统一收录在 [FAQ：深度学习性质预测的参考文献](../../faq.md#deep-learning-property-prediction-references)。

## 1) RNA-seq 覆盖与表达

- **作用：** 根据 DNA 序列预测 RNA-seq 覆盖和基因表达轨道。
- **预测方法：** 使用 AlphaGenome all-folds 官方权重直接推理 `RNA_SEQ` 输出；短序列会居中补 `N` 到官方支持的上下文长度，可用 `ontology_terms` 限定组织或细胞类型。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="rna_seq")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也可改用 `VariantContext` 和 `predict_variant_effects()`，得到 REF、ALT、ALT−REF 三组轨道。支持人或小鼠。
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

## 2) CAGE 转录起始与表达

- **作用：** 预测 CAGE 信号，用于描述加帽转录本的转录起始位置和表达强度。
- **预测方法：** 使用 AlphaGenome 官方 `CAGE` 输出头直接生成位置 × 轨道矩阵；`ontology_terms` 可筛选目标组织或细胞类型。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="cage")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext` 形式的单 SNV 上下文。序列最长 1,048,576 bp，支持人或小鼠。
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

## 3) PRO-cap 新生转录起始

- **作用：** 预测 PRO-cap 新生 RNA 5′ 端信号，用于定位活跃转录起始事件。
- **预测方法：** 使用 AlphaGenome 官方 `PROCAP` 输出头直接推理；输出分辨率记录在 `metadata["resolution_bp"]` 中。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="procap")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可指定 `organism` 和 `ontology_terms`。
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

## 4) ATAC-seq 染色质可及性

- **作用：** 预测 ATAC-seq 开放染色质轨道，用于识别可能处于可及状态的调控区域。
- **预测方法：** 使用 AlphaGenome 官方 `ATAC` 输出头直接推理，并按位置返回不同组织或细胞类型的可及性信号。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="atac")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持只包含一个碱基变化的 `VariantContext`。可通过 `ontology_terms` 限制轨道。
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

## 5) DNase-seq 染色质可及性

- **作用：** 预测 DNase-seq 超敏信号，用于描述染色质可及性。
- **预测方法：** 使用 AlphaGenome 官方 `DNASE` 输出头直接推理；短输入会居中补齐，补齐长度保存在结果元数据中。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="dnase")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可选 `organism="human"` 或 `"mouse"`。
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

## 6) 组蛋白修饰

- **作用：** 预测组蛋白 ChIP-seq 轨道，用于表征不同组蛋白修饰在序列区域中的信号。
- **预测方法：** 使用 AlphaGenome 官方 `CHIP_HISTONE` 输出头直接推理，返回位置 × 组蛋白轨道矩阵。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="chip_histone")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可按物种和 ontology term 选择输出。
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

## 7) 转录因子结合

- **作用：** 预测转录因子 ChIP-seq 结合轨道，用于定位可能的转录因子结合信号。
- **预测方法：** 使用 AlphaGenome 官方 `CHIP_TF` 输出头直接推理，返回位置 × 转录因子轨道矩阵。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="chip_tf")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可用 `ontology_terms` 限定组织或细胞类型。
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

## 8) 剪接供体/受体位点

- **作用：** 预测序列各位置的剪接供体和剪接受体相关轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_SITES` 输出头直接推理，并保留输出分辨率及轨道名称。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_sites")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`，变异结果额外包含 ALT−REF 差值轴。
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

## 9) 剪接位点使用比例

- **作用：** 预测候选剪接位点的使用信号或使用比例轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_SITE_USAGE` 输出头直接推理；具体轨道身份由 `output_names` 给出。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_site_usage")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可指定物种、ontology term 和官方上下文长度。
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

## 10) 剪接连接计数

- **作用：** 预测候选 splice junction 的连接计数轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_JUNCTIONS` 输出头直接推理；结果元数据中的 `junctions` 保存候选连接标识。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_junctions")[必须], backend[可选])`
- **输入：** 只接受 `BiologicalSequence`。该任务的候选连接集合可能随等位基因变化，因此当前不接受 `VariantContext`。
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

## 11) 三维 DNA 接触图

- **作用：** 预测染色质区域之间的三维 DNA 接触概率图。
- **预测方法：** 使用 AlphaGenome 官方 `CONTACT_MAPS` 输出头直接推理，输出两个位置轴和一个轨道轴。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="contact_maps")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`，变异结果在最前面增加 REF、ALT、ALT−REF 轴。
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

## 12) 人类调控轨道

- **作用：** 预测人类基因组的 5,313 条 CAGE、DNase/ATAC、转录因子 ChIP 和组蛋白 ChIP 轨道。
- **预测方法：** 使用公开的 Enformer 预训练权重直接推理；默认把输入居中补 `N` 到 196,608 bp，并输出 896 个 128 bp bin。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="enformer", task="human_tracks")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也可用 `VariantContext` 和 `predict_variant_effects()`，获得 REF、ALT、ALT−REF。`max_length` 必须是 114,688–393,216 范围内的 128 倍数。
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

## 13) 小鼠调控轨道

- **作用：** 预测小鼠基因组的 1,643 条 CAGE、DNase/ATAC、转录因子 ChIP 和组蛋白 ChIP 轨道。
- **预测方法：** 使用 Enformer 的小鼠输出头直接推理；默认上下文为 196,608 bp，输出为 128 bp 分箱轨道。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="enformer", task="mouse_tracks")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持只含一个 SNV 的 `VariantContext`。输入不足上下文长度时自动在两端居中补 `N`。
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

## 14) 14 类单碱基基因组分割

- **作用：** 为人 DNA 序列中的每个碱基预测 14 类基因和调控元件概率。
- **预测方法：** 使用带训练完成分割头的 `InstaDeepAI/segment_nt`，即 NT-v2 骨干网络加 U-Net 分割头；对二分类 logits 做 softmax 后返回“该特征存在”的概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="segmentnt", task="genomic_segmentation", allow_remote_code=True)[必须], backend[可选])`
- **输入：** 最长 30,000 bp 的人 `BiologicalSequence`。14 类为 protein-coding gene、lncRNA、exon、intron、splice donor、splice acceptor、5′ UTR、3′ UTR、CTCF-bound、polyA signal，以及组织特异/不变的 enhancer 和 promoter。
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

## 15) 长上下文零样本变异效应

- **作用：** 计算单 SNV 的 REF、ALT 序列似然及 ALT−REF 变异效应分数。
- **预测方法：** 使用 Evo 2 7B 官方基础模型分别计算完整 REF 和 ALT 上下文的平均序列 log-likelihood，再计算 `alternate − reference`。该值是零样本排序分数，不是临床校准的致病概率。
- **API：** `dnakit.predictions.predict_variant_effects(inputs[必须], config=PropertyPredictionConfig(model="evo2", task="variant_effect")[必须], backend[可选])`
- **输入：** `VariantContext`；REF 与 ALT 必须等长且只有一个 A/C/G/T 位点不同。
- **示例代码：**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="A" * 4096 + "C" + "G" * 4095,
    alternate_sequence="A" * 4096 + "T" + "G" * 4095,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="evo2",
        task="variant_effect",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **示例结果：**

```text
('reference_likelihood', 'alternate_likelihood', 'alternate_minus_reference')
(3,)
```

## 16) 外显子概率

- **作用：** 根据一个候选位点两侧的正向和反向上下文预测该位点属于 exon 的概率。
- **预测方法：** 从 Evo 2 `blocks.26` 读取正向、反向上下文的末 token 4,096 维 embedding，拼接为 8,192 维，再通过官方已经训练好的 MLP 头输出概率；不会重新训练分类器。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="evo2", task="exon_probability", allow_remote_code=True)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`，顺序必须是 forward context、reverse context；两条都使用 `sequence_type="gene"`。
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

## 17) 等位基因条件概率变异效应

- **作用：** 根据变异位点上游序列，比较下一 token 支持 REF 与 ALT 碱基的概率。
- **预测方法：** 使用 GENERator v2 eukaryote 1.2B 官方权重；最多读取 8,192 bp 上游上下文，去除开头 `N` 并左截断到 6 的倍数，然后计算 `log(p_ref/(p_alt+1e-10))`。该值是零样本排序分数，不是临床校准概率。
- **API：** `dnakit.predictions.predict_variant_effects(inputs[必须], config=PropertyPredictionConfig(model="generator", task="variant_effect", allow_remote_code=True)[必须], backend[可选])`
- **输入：** `VariantContext`；REF 与 ALT 必须等长且只含一个 SNV。模型只使用变异位点之前的参考上下文。
- **示例代码：**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="ACGT" * 1024 + "C" + "G" * 255,
    alternate_sequence="ACGT" * 1024 + "T" + "G" * 255,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="generator",
        task="variant_effect",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **示例结果：**

```text
('reference_probability', 'alternate_probability', 'log_ratio_score')
(3,)
```

## 18) 核酸–蛋白中心法则关系

- **作用：** 预测一条 gene 序列与一条 protein 序列是否具有 CentralDogma 关系。
- **预测方法：** 调用 LucaOneTasks 官方 `CentralDogma` 二分类 checkpoint 和 `src/predict_v1.py` 直接推理，返回正类概率及官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="central_dogma", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；第一条必须为 `sequence_type="gene"`，第二条必须为 `sequence_type="protein"`。`model_source_path` 必须指向 LucaOneTasks 官方源码目录。
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

## 19) SupKTax 分类

- **作用：** 使用官方 SupKTax checkpoint 对 gene 序列进行分类。
- **预测方法：** 调用 LucaOneTasks 官方 `SupKTax` 多分类 checkpoint 直接推理，并返回概率最高的 `top_k` 个 `label.txt` 标签。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="supktax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`。结果标签以该 checkpoint 的 `label.txt` 为准。
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

## 20) GenusTax 分类

- **作用：** 使用官方 GenusTax checkpoint 对 gene 序列进行分类。
- **预测方法：** 调用 LucaOneTasks 官方 `GenusTax` 多分类 checkpoint 直接推理，并按概率返回 Top-k 标签。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="genustax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`。LucaOneTasks README 对 SupKTax/GenusTax 的自然语言层级存在对调，因此解释结果时以实际 `label.txt` 为准。
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

## 21) 物种分类

- **作用：** 预测 gene 序列对应的物种类别。
- **预测方法：** 调用 LucaOneTasks 官方 `SpeciesTax` 多分类 checkpoint 直接推理，返回 Top-k 物种标签和对应概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="speciestax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`；可用 `top_k` 设置返回标签数。
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

## 22) 原核蛋白亚细胞定位

- **作用：** 预测原核蛋白的亚细胞定位类别。
- **预测方法：** 调用 LucaOneTasks 官方 `ProtLoc` 多分类 checkpoint 直接推理，返回 Top-k 定位标签及概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_location", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="protein"` 的 `BiologicalSequence`。
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

## 23) 蛋白稳定性

- **作用：** 根据蛋白序列输出官方 ProtStab 任务定义下的稳定性回归值。
- **预测方法：** 调用 LucaOneTasks 官方 `ProtStab` 回归 checkpoint 和 `predict_v1.py` 直接推理，不在 DNAKit 中拟合回归器。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_stability", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="protein"` 的 `BiologicalSequence`。
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
        task="protein_stability",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **示例结果（数值取决于输入和 checkpoint）：**

```text
('prediction',)
(1,)
```

## 24) ncRNA 家族

- **作用：** 预测非编码 RNA 序列所属的家族。
- **预测方法：** 调用 LucaOneTasks 官方 `ncRNAFam` 多分类 checkpoint 直接推理，返回 Top-k 家族标签和概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="ncrna_family", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`；DNA 或 RNA 字母均可，`U` 由官方任务流程处理。
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

## 25) 甲型流感抗原关系

- **作用：** 预测两条甲型流感 gene 序列之间的抗原关系类别。
- **预测方法：** 调用 LucaOneTasks 官方 `InfA` 二分类 checkpoint 直接推理，返回正类概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="influenza_antigenicity", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；两条序列都必须是 `sequence_type="gene"`。
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

## 26) 蛋白–蛋白互作

- **作用：** 预测两条蛋白序列是否具有蛋白–蛋白相互作用。
- **预测方法：** 调用 LucaOneTasks 官方 `PPI` 二分类 checkpoint 直接推理，返回互作概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_interaction", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；两条序列都必须是 `sequence_type="protein"`。
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

## 27) ncRNA–蛋白互作

- **作用：** 预测一条 ncRNA 与一条蛋白序列是否具有相互作用。
- **预测方法：** 调用 LucaOneTasks 官方 `ncRPI` 二分类 checkpoint 直接推理，返回互作概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="ncrna_protein_interaction", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；第一条必须是 `sequence_type="gene"`，第二条必须是 `sequence_type="protein"`。
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
