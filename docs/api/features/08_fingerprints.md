# 序列表征与特征处理

将 DNA 序列转为数值向量，进行特征缩放、降维与筛选。

## 1) 整数编码 {#1-fp-001}

- **作用：** 按显式码表把每个 DNA 字符转换成整数，同时保留原始位置顺序，便于分词、模型输入和可逆检查编码规则。
- **计算方法：** 默认按 A=0、C=1、G=2、T=3 逐字符编码，其他符号按显式策略处理。
- **API：** `dnakit.fingerprints.integer_encode(value[必须], ambiguity_policy[可选], gap_policy[可选], max_output_length[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 IUPAC 策略、Gap 策略和最大输出长度。
- **输出：** 返回 `IntegerEncodingResult`，示例展示 `values`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import integer_encode

result = integer_encode(DNASequence("ACGT"))
print(result.values)
```

- **示例结果：**

```text
(0, 1, 2, 3)
```

## 2) One-hot编码 {#2-fp-002-one-hot}

- **作用：** 把每个碱基转换为对应 A/C/G/T 的四维向量，形成“序列长度 × 4”的位置保持矩阵，便于统计模型和神经网络使用。
- **计算方法：** 按 base_order 生成每个位置的四列向量；标准碱基对应列为 1，其余列为 0。
- **API：** `dnakit.fingerprints.one_hot_encode(value[必须], ambiguity_policy[可选], gap_policy[可选], base_order[可选], max_output_length[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 IUPAC、Gap、列顺序和输出长度策略。
- **输出：** 返回 `OneHotEncodingResult`，示例展示 `values`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import one_hot_encode

result = one_hot_encode(DNASequence("ACGT"))
print(result.values[0], result.values[-1])
```

- **示例结果：**

```text
(1.0, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 1.0)
```

## 3) k-mer 特征 {#3-fp-003-k-mer}

- **作用：** 把全部 k-mer 的数量、频率或存在性按固定顺序组成特征向量，用于序列比较和传统机器学习；通过 `canonical=True` 可把反向互补 k-mer 合并到同一特征。
- **计算方法：** 滑窗统计 k-mer；count 保留计数，frequency 归一化，binary 仅记录存在性。
- **API：** 推荐使用 `dnakit.fingerprints.kmer(value[必须], k[必须], canonical[可选], mode[可选], representation[可选], ambiguity_policy[可选], overlapping[可选], cross_gaps[可选], max_dimension[可选])`；旧名称 `kmer_fingerprint(...)` 仅保留为兼容别名。
- **输入：** 必填序列和 `k`；可选 `canonical`、`mode="count"|"frequency"|"binary"`、存储表示、重叠策略、`ambiguity_policy="error"|"ignore"` 和 `cross_gaps`。
- **输出：** 返回 k-mer 特征结果，包含特征值 `values` 和维度 `dimension`。
- **术语：** `count`、`frequency` 和 `binary` 均归入 k-mer 特征；`binary` 是存在性特征。需要固定长度 DNA 指纹时，请使用 [Hashed k-mer 位指纹](08_feature_engineering.md#1-hashed-k-mer)。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import kmer

sequence = DNASequence("ACGT")
standard = kmer(sequence, k=2, representation="sparse")
canonical = kmer(sequence, k=2, canonical=True, representation="sparse")
print(dict(standard.values), standard.dimension)
print(dict(canonical.values), canonical.dimension)
```

- **示例结果：**

```text
{'AC': 1, 'CG': 1, 'GT': 1} 16
{'AC': 2, 'CG': 1} 10
```

## 4) k-mer Sketch（MinHash/FracMinHash） {#4-fp-005-k-mer-sketchminhashfracminhash}

- **作用：** 对序列的 k-mer 哈希进行 Bottom-k 或阈值抽样，生成体积较小且可复现的 Sketch，用于近似估计 Jaccard 相似度和检索大批序列。
- **计算方法：** MinHash 保留最小的 num_hashes 个唯一哈希；FracMinHash 按 scaled 阈值保留哈希。
- **API：** `dnakit.fingerprints.minhash(value[必须], k[必须], num_hashes[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`、`dnakit.fingerprints.fracminhash(value[必须], k[必须], scaled[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`
- **输入：** 必填序列和 `k`；MinHash 可选 `num_hashes`，FracMinHash 可选 `scaled`，两者可选 seed/canonical/上限。
- **输出：** 返回 `SketchResult`，包含哈希值和 sketch 参数。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import minhash

sequence = DNASequence("ACGTACGT")
first = minhash(sequence, k=3, num_hashes=3, seed=7)
second = minhash(sequence, k=3, num_hashes=3, seed=7)
print(first.hashes == second.hashes, len(first.hashes))
```

- **示例结果：**

```text
True 2
```

<a id="neural-representations"></a>

## 5) 神经网络表征

- **作用：** 将每条 DNA 序列转换为固定长度浮点向量，供聚类、检索或机器学习使用。
- **计算方法：** 按 mean/cls/max/last 汇聚模型的 token 表征；超过上下文长度时分块，再对块向量等权平均。
- **API：** `dnakit.representations.extract_representations(records, config=RepresentationConfig(...))`。
- **输入：** 非空、无显式 Gap 的 DNA 记录；可选模型、checkpoint、汇聚方式和设备。默认 LucaOne 后端需要模型依赖及 `allow_remote_code=True`。
- **输出：** 返回 `RepresentationResult`，示例展示 `model_name`、`representations`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.representations import RepresentationConfig, extract_representations

records = [
    DNARecord(DNASequence("ACGTACGT"), "seq-1"),
    DNARecord(DNASequence("AACCGGTT"), "seq-2"),
]
result = extract_representations(
    records,
    config=RepresentationConfig(
        model="lucaone",
        pooling="mean",
        allow_remote_code=True,
    ),
)
print(result.model_name)
print(result.representations.shape)
```

- **示例结果：** `model_name` 为 `lucaone`；矩阵含 2 行，列数取决于 checkpoint。本示例未在此处加载模型权重实测。

[FAQ：补充说明](../../faq/tutorial_basics.md#embedding-models) · [FAQ：补充说明](../../faq/tutorial_basics.md#embedding-loading)

<span id="51-api"></span>
<span id="52-11-rep"></span>
<span id="53-rep"></span>

## 6) 把 DNA 序列编码成特征矩阵 {#bio-23}

- **作用：** 得到可用于分类、聚类或统计的数值特征。
- **计算方法：** 按指定 descriptor 计算；示例 NAC 计算 A/C/G/T 频率，RCKmer 合并反向互补词频。
- **API：** `dnakit.bioinformatics.encoding.encode_dna(sequences, descriptor, parameters=None, property_tables=None, labels=None, training_ids=None)`
- **输入：** `DNA([...])` 或 ID→A/C/G/T 序列字典；需要 `bio-analysis`。位置编码要求等长，理化编码需显式性质表；[全部编码与参数](../../faq/bioinformatics.md#encoding)。
- **输出：** 返回 `EncodingResult`，包含序列 ID `ids` 和数值特征矩阵 `values`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.bioinformatics.encoding import encode_dna

sequences = DNA(["AACGTG", "CCGTAT"])
nac = encode_dna(sequences, "NAC")
rckmer = encode_dna(sequences, "RCKmer", parameters={"k": 2})
print(nac.ids)
print(nac.values.round(3).tolist())
print(rckmer.values.shape)
```

- **示例结果：**

```text
('sequence_1', 'sequence_2')
[[0.333, 0.167, 0.333, 0.167], [0.167, 0.333, 0.167, 0.333]]
(2, 10)
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#encoding)。

## 7) 特征缩放、降维与筛选 {#bio-24}

- **作用：** 拟合训练集特征转换，并对新样本复用。
- **计算方法：** 按 `scale` 缩放后运行 PCA 等方法；监督筛选用 `labels`，新样本调用 `transform`。
- **API：** `dnakit.bioinformatics.features.analyze_features(values, method='pca', labels=None, feature_names=None, scale='standard', components=2, selected_features=10)`
- **输入：** samples×features 矩阵；需要 `bio-analysis`。训练与新样本列顺序一致；t-SNE 无 `transform`；[方法清单](../../faq/bioinformatics.md#encoding)。
- **输出：** 返回 `FeatureAnalysis`，包含转换后的 `values`；`transform` 可处理新样本。
- **示例代码：**

```python
from dnakit.bioinformatics.features import analyze_features

fitted = analyze_features([[1, 2], [2, 1], [3, 4]],
                          feature_names=["x", "y"], method="pca", components=1)
held_out = fitted.transform([[2, 3]], feature_names=["x", "y"])
print(fitted.values.shape, held_out.shape)
```

- **示例结果：**

```text
(3, 1) (1, 1)
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#encoding)。
