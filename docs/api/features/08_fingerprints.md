# DNA表征

将 DNA 序列embedding得到表征。

## 1) FP-001 整数编码

- **作用：** 按显式码表把每个 DNA 字符转换成整数，同时保留原始位置顺序，便于分词、模型输入和可逆检查编码规则。
- **API：** `dnakit.fingerprints.integer_encode(value[必须], ambiguity_policy[可选], gap_policy[可选], max_output_length[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 IUPAC 策略、Gap 策略和最大输出长度。
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

## 2) FP-002 One-hot编码

- **作用：** 把每个碱基转换为对应 A/C/G/T 的四维向量，形成“序列长度 × 4”的位置保持矩阵，便于统计模型和神经网络使用。
- **API：** `dnakit.fingerprints.one_hot_encode(value[必须], ambiguity_policy[可选], gap_policy[可选], base_order[可选], max_output_length[可选])`
- **输入：** 必填 `DNASequence`/`DNARecord`；可选 IUPAC、Gap、列顺序和输出长度策略。
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

## 3) FP-003 k-mer 特征

- **作用：** 把全部 k-mer 的数量、频率或存在性按固定顺序组成特征向量，用于序列比较和传统机器学习；通过 `canonical=True` 可把反向互补 k-mer 合并到同一特征。
- **API：** 推荐使用 `dnakit.fingerprints.kmer(value[必须], k[必须], canonical[可选], mode[可选], representation[可选], ambiguity_policy[可选], overlapping[可选], cross_gaps[可选], max_dimension[可选])`；旧名称 `kmer_fingerprint(...)` 仅保留为兼容别名。
- **输入：** 必填序列和 `k`；可选 `canonical`、`mode="count"|"frequency"|"binary"`、存储表示、重叠策略、`ambiguity_policy="error"|"ignore"` 和 `cross_gaps`。
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

## 4) FP-005 k-mer Sketch（MinHash/FracMinHash）

- **作用：** 对序列的 k-mer 哈希进行 Bottom-k 或阈值抽样，生成体积较小且可复现的 Sketch，用于近似估计 Jaccard 相似度和检索大批序列。
- **API：** `dnakit.fingerprints.minhash(value[必须], k[必须], num_hashes[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`、`dnakit.fingerprints.fracminhash(value[必须], k[必须], scaled[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`
- **输入：** 必填序列和 `k`；MinHash 可选 `num_hashes`，FracMinHash 可选 `scaled`，两者可选 seed/canonical/上限。
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

## 5) 神经网络表征

- **作用：** 使用预训练 DNA 基础模型把每条序列压缩为固定长度浮点向量，用于聚类、检索或下游机器学习；向量含义取决于模型、checkpoint 和 pooling 配置。

神经网络表征是 `DATA-027` 的 rep 提取部分：使用预训练 DNA 基础模型，将每条 DNA 序列转换为一个固定长度的 float32 向量。将 rep 用于 k-means 的流程见 [DATA-027 神经网络聚类](10_clustering.md#data-027-neural-clustering)。

### 5.1 表征 API

- **API：** `dnakit.representations.extract_representations(records[必须], config[可选])`。
- **配置：** `dnakit.representations.RepresentationConfig`，可选模型、checkpoint、pooling、设备、dtype、批大小、最大长度和进度条。
- **输入：** 一条或多条 `DNARecord`；不接受空序列或显式 `Gap`。
- **输出：** `RepresentationResult`，其 `representations` 是只读二维矩阵，每行对应一条输入记录；实际维度保存在 `embedding_dimension`。
- **汇聚规则：** 支持 `mean`、`cls`、`max`、`last` pooling。超过模型上下文窗口的序列会分块，各块池化向量最后等权平均。

### 5.2 支持的 11 种神经网络的 rep

以下 `model` 名称直接来自当前 `MODEL_REGISTRY`：

| 序号 | `model`           | rep 来源                                     | 默认官方 checkpoint                                                                                                                        | 加载条件                                                                                                   |
| ---: | ------------------- | -------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------- |
|    1 | `alphagenome`     | AlphaGenome all folds                        | [google/alphagenome-all-folds](https://huggingface.co/google/alphagenome-all-folds)                                                         | 须接受模型条款并安装[官方 research 代码](https://github.com/google-deepmind/alphagenome_research)           |
|    2 | `caduceus`        | Caduceus-Ph 131k                             | [kuleshov-group/caduceus-ph 131k](https://huggingface.co/kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16)                     | `neural`、`neural-caduceus` extras；须显式 `allow_remote_code=True`                                  |
|    3 | `dnabert2`        | DNABERT-2-117M                               | [zhihan1996/DNABERT-2-117M](https://huggingface.co/zhihan1996/DNABERT-2-117M)                                                               | `neural` extra；须显式 `allow_remote_code=True`                                                        |
|    4 | `enformer`        | Enformer PyTorch                             | [EleutherAI/enformer-official-rough](https://huggingface.co/EleutherAI/enformer-official-rough)                                             | `neural`、`neural-enformer` extras                                                                     |
|    5 | `evo2`            | Evo 2 7B                                     | [arcinstitute/evo2_7b](https://huggingface.co/arcinstitute/evo2_7b)                                                                         | `neural`、`neural-evo2` extras；需兼容 GPU 环境                                                        |
|    6 | `generator`       | GENERator v2 eukaryote 1.2B                  | [GenerTeam/GENERator-v2-eukaryote-1.2b-base](https://huggingface.co/GenerTeam/GENERator-v2-eukaryote-1.2b-base)                             | `neural` extra；须显式 `allow_remote_code=True`                                                        |
|    7 | `grover`          | GROVER                                       | [PoetschLab/GROVER](https://huggingface.co/PoetschLab/GROVER)                                                                               | `neural` extra；已完成真实 checkpoint rep + k-means smoke                                                |
|    8 | `hyenadna`        | HyenaDNA medium 450k                         | [LongSafari/hyenadna-medium-450k-seqlen-hf](https://huggingface.co/LongSafari/hyenadna-medium-450k-seqlen-hf)                               | `neural` extra；须显式 `allow_remote_code=True`                                                        |
|    9 | `janusdna`        | JanusDNA 131k                                | [Harvard Dataverse DOI 10.7910/DVN/HDT0RN](https://doi.org/10.7910/DVN/HDT0RN)                                                              | checkpoint 校验官方 MD5；还须提供[官方源码](https://github.com/Qihao-Duan/JanusDNA)及 `model_source_path` |
|   10 | `lucaone`（默认） | LucaOne gene step 36.8M                      | [LucaGroup/LucaOne-gene-step36.8M](https://huggingface.co/LucaGroup/LucaOne-gene-step36.8M)                                                 | `neural` extra；须显式 `allow_remote_code=True`                                                        |
|   11 | `ntv2`            | Nucleotide Transformer v2 500M multi-species | [InstaDeepAI/nucleotide-transformer-v2-500m-multi-species](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species) | `neural` extra；须显式 `allow_remote_code=True`                                                        |

查询当前注册的全部模型，不会下载 checkpoint：

```python
from dnakit.representations import available_embedding_models

print(available_embedding_models())
```

```text
('alphagenome', 'caduceus', 'dnabert2', 'enformer', 'evo2', 'generator', 'grover', 'hyenadna', 'janusdna', 'lucaone', 'ntv2')
```

### 5.3 提取 rep 示例

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

默认 checkpoint 下载到当前工作目录的 `ckpt/lucaone-gene-step36-8m/`，已有完整文件时直接复用。下载和逐序列提取默认显示进度条；可用 `checkpoint_dir`、`checkpoint_path` 或 `show_progress` 显式调整。
