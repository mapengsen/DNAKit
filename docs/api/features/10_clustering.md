# 聚类

根据序列一致性、k-mer、指纹相似度或 DNA 基础模型 rep 对序列进行聚类，从而得到序列分组及代表序列。

## 1) DATA-007 · Identity聚类

**作用：** 计算序列两两相似度并按阈值图的连通分量分组，返回 cluster 标签、成员和代表记录，用于整理相近序列及观察数据结构。


**API：** `dnakit.datasets.cluster_sequences(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.ClusterConfig`，本项设置 `method="identity"`。

**输入：** 必填记录集合；可选 identity 阈值、代表策略和比对资源上限。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, cluster_sequences

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = cluster_sequences(
    records, config=ClusterConfig(method="identity", threshold=0.7)
)
print(result.labels)
```

**示例结果：**

```text
(0, 0, 1)
```

**限制：** 采用有界 exact pairwise alignment；阈值图的连通关系不等于所有组内两两均过阈值。

## 2) DATA-008 · k-mer聚类

**作用：** 根据序列之间的 k-mer Jaccard、Containment 或 Cosine 相似度建立分组，返回 cluster 成员和标签，适合不做碱基级比对的组成型聚类。


**API：** `dnakit.datasets.cluster_sequences(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.ClusterConfig`，本项设置 `method="kmer"`。

**输入：** 必填记录集合；可选 k、canonical、阈值、代表策略和两两比较上限。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, cluster_sequences

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = cluster_sequences(
    records, config=ClusterConfig(method="kmer", threshold=0.5, k=2)
)
print(result.labels)
```

**示例结果：**

```text
(0, 0, 1)
```

**限制：** 这是原生穷举阈值聚类，不自动调用 Mash、sourmash 或 Dashing。

## 3) DATA-009 · 指纹聚类

**作用：** 使用固定 schema DNA 指纹的 Tanimoto/Jaccard 等相似度对序列分组，返回 cluster 标签和成员，用于按局部模式或 Panel 特征聚类。


**API：** `dnakit.datasets.cluster_sequences(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.ClusterConfig`，本项设置 `method="fingerprint"`。

**输入：** 必填记录集合；可选阈值、代表策略和资源上限。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, cluster_sequences

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = cluster_sequences(
    records, config=ClusterConfig(method="fingerprint", threshold=0.5)
)
print(result.labels, result.method)
```

**示例结果：**

```text
(0, 0, 1) fingerprint
```

**限制：** 当前 API 使用内置确定性指纹配置，不接收任意外部预计算矩阵。

## 4) DATA-010 · 层次聚类

**作用：** 根据预先计算的序列距离逐步合并最近簇，返回 linkage、层次关系和切分标签，用于绘制树状图及观察不同距离尺度的分组。


**API：** `dnakit.datasets.hierarchical_cluster(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.HierarchicalClusteringConfig`。

**输入：** 必填记录集合；可选 identity/edit/k-mer/fingerprint 方法和 single/complete/average linkage。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import HierarchicalClusteringConfig, hierarchical_cluster

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = hierarchical_cluster(
    records, config=HierarchicalClusteringConfig(linkage="average")
)
print(len(result.linkage), result.linkage[-1].member_count)
```

**示例结果：**

```text
2 3
```

**限制：** 返回 linkage 数据而非绘图；有界全矩阵实现适合中小数据集。

## 5) DATA-011 · 代表序列选择

**作用：** 按首条、最长、medoid 等明确策略从每个 cluster 选择一条代表序列，并保留成员映射，用于压缩数据集和后续人工检查。


**API：** `dnakit.datasets.select_representatives(records[必须], labels[必须], policy[可选], medoid_method[可选], k[可选], canonical[可选], max_records[可选], max_pairwise_comparisons[可选], max_alignment_cells[可选])`。

**输入：** 必填记录和等长 labels；可选选择策略、medoid 方法和资源上限。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import select_representatives

records = [
    DNARecord(DNASequence("AAAA"), "long"),
    DNARecord(DNASequence("AA"), "short"),
    DNARecord(DNASequence("CCCC"), "other"),
]
result = select_representatives(records, [7, 7, 2], policy="shortest")
print(result.representative_ids)
```

**示例结果：**

```text
('other', 'short')
```

**限制：** best-quality 依赖 `phred_quality`；medoid 会触发有界两两计算。

## 6) DATA-027 · 神经网络聚类 {#data-027-neural-clustering}

**作用：** 使用所选 DNA 基础模型为每条序列提取固定长度 rep，再对向量执行归一化、可选 PCA 和固定 seed 的 k-means，返回标签、中心、评分及中心最近代表序列，用于按模型表征分组。

这里使用的是 **k-means**，不是 k-mer 相似度聚类。若需要不加载神经网络的
k-mer 聚类，请使用本页的 `DATA-008`。

**API：**

- `dnakit.representations.extract_representations(records[必须], config[可选])`：只提取 rep；
- `dnakit.datasets.neural_cluster_sequences(records[必须], config[可选])`：提取 rep 后执行 k-means；
- 表征配置使用 `dnakit.representations.RepresentationConfig`；聚类配置使用
  `dnakit.datasets.NeuralClusteringConfig`。

可独立使用的 rep API、全部 11 种模型及其 checkpoint/依赖边界另见
[序列表征页中的神经网络表征栏目](08_fingerprints.md#neural-representations)。

**输入：** 必填 `DNARecord` 集合；可选模型、checkpoint 位置、mean/cls/max/last
池化、设备、dtype、批大小、L2 归一化、PCA 维数、聚类数和随机种子。显式
`Gap` 不支持；IUPAC 模糊碱基默认映射为 `N`，也可设置为直接报错。

**聚类流程：** 每条序列按模型上下文窗口分块，先得到各块的池化向量，再对块
向量求平均；随后默认做 L2 归一化，可选 PCA，最后执行固定 seed 的
`k-means++`/Lloyd 聚类。每组代表序列取离聚类中心最近的输入记录；中心位于
归一化及可选 PCA 后的空间。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import NeuralClusteringConfig, neural_cluster_sequences
from dnakit.representations import RepresentationConfig

records = [
    DNARecord(DNASequence("AAAAAA"), "a"),
    DNARecord(DNASequence("AAAAAT"), "b"),
    DNARecord(DNASequence("CCCCCC"), "c"),
    DNARecord(DNASequence("CCCCCG"), "d"),
]
result = neural_cluster_sequences(
    records,
    config=NeuralClusteringConfig(
        representation=RepresentationConfig(
            batch_size=2,
            allow_remote_code=True,
        ),
        n_clusters=2,
        seed=7,
    ),
)
print(result.labels)
print(result.representatives.ids)
print(result.checkpoint_path)
```

默认模型是 `lucaone`，对应 checkpoint
`LucaGroup/LucaOne-gene-step36.8M`。首次调用默认下载到运行命令所在目录的
`ckpt/lucaone-gene-step36-8m/`；完整 checkpoint 已存在时直接复用，不重复下载。
每个目录会保存 `.dnakit-checkpoint.json` 来源清单。可用 `checkpoint_dir` 修改缓存
根目录，或用 `checkpoint_path` 指向现有 checkpoint 目录并完全跳过下载。下载和
逐序列提取默认显示进度条。

LucaOne checkpoint 包含自定义 Transformers 代码。DNAKit 不会隐式执行它；使用
默认模型的标准后端时，仍须在 `RepresentationConfig` 中显式设置
`allow_remote_code=True`。

只提取 rep：

```python
from dnakit.representations import RepresentationConfig, extract_representations

reps = extract_representations(
    records,
    config=RepresentationConfig(allow_remote_code=True),
)
print(reps.representations.shape)
```

<span id="checkpoint"></span>**可选模型与官方 checkpoint**

| `model` | 默认官方 checkpoint | 加载条件 |
| --- | --- | --- |
| `grover` | [PoetschLab/GROVER](https://huggingface.co/PoetschLab/GROVER) | `neural` extra；本地独立模型环境已完成真实 checkpoint 的 rep + k-means smoke |
| `dnabert2` | [zhihan1996/DNABERT-2-117M](https://huggingface.co/zhihan1996/DNABERT-2-117M) | `neural` extra；需显式 `allow_remote_code=True` |
| `ntv2` | [InstaDeepAI/nucleotide-transformer-v2-500m-multi-species](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species) | `neural` extra；需显式 `allow_remote_code=True` |
| `hyenadna` | [LongSafari/hyenadna-medium-450k-seqlen-hf](https://huggingface.co/LongSafari/hyenadna-medium-450k-seqlen-hf) | `neural` extra；需显式 `allow_remote_code=True` |
| `caduceus` | [kuleshov-group/caduceus-ph 131k](https://huggingface.co/kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16) | `neural`、`neural-caduceus` extras；需显式 `allow_remote_code=True` |
| `lucaone`（默认） | [LucaGroup/LucaOne-gene-step36.8M](https://huggingface.co/LucaGroup/LucaOne-gene-step36.8M) | `neural` extra；需显式 `allow_remote_code=True` |
| `generator` | [GenerTeam/GENERator-v2-eukaryote-1.2b-base](https://huggingface.co/GenerTeam/GENERator-v2-eukaryote-1.2b-base) | `neural` extra；需显式 `allow_remote_code=True` |
| `enformer` | [EleutherAI/enformer-official-rough](https://huggingface.co/EleutherAI/enformer-official-rough) | `neural`、`neural-enformer` extras |
| `alphagenome` | [google/alphagenome-all-folds](https://huggingface.co/google/alphagenome-all-folds) | 须接受模型条款并安装[官方 research 代码](https://github.com/google-deepmind/alphagenome_research)；硬件和 JAX/Orbax 环境要求高 |
| `janusdna` | [Harvard Dataverse DOI 10.7910/DVN/HDT0RN](https://doi.org/10.7910/DVN/HDT0RN) | checkpoint 自动校验官方 MD5；还须取得[官方源码](https://github.com/Qihao-Duan/JanusDNA)并设置 `model_source_path` |
| `evo2` | [arcinstitute/evo2_7b](https://huggingface.co/arcinstitute/evo2_7b) | `neural`、`neural-evo2` extras；遵循[官方 Evo 2 环境](https://github.com/ArcInstitute/evo2)和 GPU 要求 |

除 GROVER 真实 smoke 外，其余模型在本功能中已完成注册、checkpoint 解析及对应
adapter 接入，但尚未在同一台机器上逐个完成全 checkpoint 数值验证。因此
`DATA-027` 状态为 `conditional`，不能把“已接入”解释为所有模型、硬件和版本均
已验证。

**远程代码边界：** 表中注明的模型会从 checkpoint 加载 Python 代码。DNAKit
默认拒绝这类执行，只有调用方审查官方仓库和 checkpoint 后显式设置
`allow_remote_code=True` 才加载。checkpoint 不随 wheel/sdist 分发；模型许可、
访问条款、显存和存储需求均由使用者另行确认。

**结果：** `NeuralClusteringResult` 保存 labels、cluster 成员、中心最近代表序列、
模型名、checkpoint 路径、原始/聚类维数、inertia、可计算时的 silhouette、PCA
解释方差、seed 和实际迭代数。

**限制：** 该功能只使用预训练模型提取 rep 并做无监督聚类，不提供启动子活性、
表达量、结合强度等任务型预测；超长序列的多个块采用等权平均。不同模型、
checkpoint、pooling、依赖版本和硬件精度得到的 rep 不可混作同一特征空间。
固定 seed 用于稳定初始化；不同 BLAS、PyTorch、设备或 dtype 下的浮点末位不保证
逐位相同。
