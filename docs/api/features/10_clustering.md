# 聚类

根据序列一致性、k-mer、指纹相似度或 DNA 基础模型 rep 对序列进行聚类，从而得到序列分组及代表序列。

## 1) Identity聚类 {#1-data-007-identity}

- **作用：** 计算序列两两相似度并按阈值图的连通分量分组，返回 cluster 标签、成员和代表记录，用于整理相近序列及观察数据结构。
- **计算方法：** 两两计算 identity，超过阈值的记录连边；取阈值图的连通分量作为簇。


- **API：** `dnakit.datasets.cluster_sequences(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.ClusterConfig`，本项设置 `method="identity"`。

- **输入：** 必填记录集合；可选 identity 阈值、代表策略和比对资源上限。
- **输出：** 返回 `ClusteringResult`，包含分组 `clusters`、标签 `labels` 和代表记录。

- **示例代码：**

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

- **示例结果：**

```text
(0, 0, 1)
```

## 2) k-mer聚类 {#2-data-008-k-mer}

- **作用：** 根据序列之间的 k-mer Jaccard、Containment 或 Cosine 相似度建立分组，返回 cluster 成员和标签，适合不做碱基级比对的组成型聚类。
- **计算方法：** 计算 k-mer 相似度并建立阈值图，按连通分量分组。


- **API：** `dnakit.datasets.cluster_sequences(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.ClusterConfig`，本项设置 `method="kmer"`。

- **输入：** 必填记录集合；可选 k、canonical、阈值、代表策略和两两比较上限。
- **输出：** 返回 `ClusteringResult`，包含分组 `clusters`、标签 `labels` 和代表记录。

- **示例代码：**

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

- **示例结果：**

```text
(0, 0, 1)
```

## 3) 指纹聚类 {#3-data-009}

- **作用：** 使用固定 schema DNA 指纹的 Tanimoto/Jaccard 等相似度对序列分组，返回 cluster 标签和成员，用于按局部模式或 Panel 特征聚类。
- **计算方法：** 计算同构指纹的相似度并建立阈值图，按连通分量分组。


- **API：** `dnakit.datasets.cluster_sequences(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.ClusterConfig`，本项设置 `method="fingerprint"`。

- **输入：** 必填记录集合；可选阈值、代表策略和资源上限。
- **输出：** 返回 `ClusteringResult`，包含分组 `clusters`、标签 `labels` 和代表记录。

- **示例代码：**

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

- **示例结果：**

```text
(0, 0, 1) fingerprint
```

## 4) 层次聚类 {#4-data-010}

- **作用：** 根据预先计算的序列距离逐步合并最近簇，返回 linkage、层次关系和切分标签，用于绘制树状图及观察不同距离尺度的分组。
- **计算方法：** 从单条记录起按 single/complete/average linkage 依次合并最近簇。


- **API：** `dnakit.datasets.hierarchical_cluster(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.HierarchicalClusteringConfig`。

- **输入：** 必填记录集合；可选 identity/edit/k-mer/fingerprint 方法和 single/complete/average linkage。
- **输出：** 返回 `HierarchicalClusteringResult`，示例展示 `linkage`。

- **示例代码：**

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

- **示例结果：**

```text
2 3
```

## 5) 代表序列选择 {#5-data-011}

- **作用：** 按首条、最长、medoid 等明确策略从每个 cluster 选择一条代表序列，并保留成员映射，用于压缩数据集和后续人工检查。
- **计算方法：** 按已有 labels 分组，再按首条、长度或 medoid 距离策略选择每组代表。


- **API：** `dnakit.datasets.select_representatives(records[必须], labels[必须], policy[可选], medoid_method[可选], k[可选], canonical[可选], max_records[可选], max_pairwise_comparisons[可选], max_alignment_cells[可选])`。

- **输入：** 必填记录和等长 labels；可选选择策略、medoid 方法和资源上限。
- **输出：** 返回 `RepresentativeSelectionResult`，示例展示 `representative_ids`。

- **示例代码：**

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

- **示例结果：**

```text
('other', 'short')
```

## 6) 神经网络聚类 {#data-027-neural-clustering}

- **作用：** 使用所选 DNA 基础模型为每条序列提取固定长度 rep，再对向量执行归一化、可选 PCA 和固定 seed 的 k-means，返回标签、中心、评分及中心最近代表序列，用于按模型表征分组。
- **计算方法：** 提取并汇聚模型向量，默认 L2 归一化，可选 PCA，再用固定 seed 的 k-means++/Lloyd 聚类；每组取距中心最近的记录作代表。
- **API：** `dnakit.datasets.neural_cluster_sequences(records, config=NeuralClusteringConfig(...))`。
- **输入：** 非空、无显式 Gap 的 DNA 记录，聚类数、seed 和模型配置。默认 LucaOne 需要对应依赖及 `allow_remote_code=True`；IUPAC 默认映射为 `N`。
- **输出：** 返回 `NeuralClusteringResult`，示例展示 `labels`、`representatives`、`checkpoint_path`。
- **示例代码：**

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

- **示例结果：** 4 个聚类标签、2 个代表序列 ID 和 checkpoint 路径；具体标签取决于模型和环境，此处未执行模型推理。

<span id="checkpoint"></span>[FAQ：补充说明](../../faq/tutorial_basics.md#neural-clustering)

## 7) 对样本聚类 {#bio-25}

- **作用：** 根据特征把相似样本分组。
- **计算方法：** 调用选定聚类方法；MCL 则对非负相似度矩阵迭代 expansion/inflation。
- **API：** `dnakit.bioinformatics.features.cluster_features(values, method='kmeans', seed=0, parameters=None)`
- **输入：** 数值特征矩阵，聚类前自行缩放；MCL 要求对称非负 affinity；需要 `bio-analysis`；[方法清单](../../faq/bioinformatics.md#encoding)。
- **输出：** 返回聚类模型；`labels_` 是每个样本的聚类标签。
- **示例代码：**

```python
from dnakit.bioinformatics.features import cluster_features

model = cluster_features([[0, 0], [0, 1], [9, 9], [9, 10]],
                         parameters={"n_clusters": 2}, seed=0)
print(model.labels_[0] == model.labels_[1])
print(model.labels_[0] != model.labels_[2])
```

- **示例结果：**

```text
True
True
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#encoding)。

## 8) 按标签进行多样性采样 {#diversity-sampling}

`dnakit.datasets.diversity_sample` 从候选池中选出指定数量的原始序列。默认使用冻结的
DNABERT2 提取 mean embedding，做 L2 归一化，再在每个类别内分别运行固定 seed 的
k-means++/Lloyd。模型只做推理，不训练。`diversity` 参数选择高或低多样性模式：

| 参数 | 聚类数量 | 选择方法 |
| --- | --- | --- |
| `diversity="high"`（默认） | 合计 `n_samples` 个簇，按标签比例分配 | 每簇选距中心最近的一条记录 |
| `diversity="low"` | 每个有采样名额的标签默认 5 簇 | 比较各中心附近的定额候选组，选最紧密的一组 |

`n_samples` 始终控制**最终采样条数**，仅在 `high` 模式下同时等于总簇数。
`low_n_clusters` 控制 `low` 模式下每个标签的簇数，默认 5。
显式传入 `n_samples`、`diversity` 或 `low_n_clusters` 会覆盖配置中的对应值，
但不修改配置对象；省略时沿用配置值，分别默认为 400、`"high"`、5。

名额按候选池的类别比例用**最大余数法**分配：先向下取整，再按余数从大到小补齐。
例如，840 条候选中标签 0 有 173 条、标签 1 有 667 条，选取 400 条时，分别分配
82 和 318 条。`high` 模式分别聚成 82 和 318 个簇，每簇选择一条；`low` 模式分别聚
5 个簇，再分别选取 82 和 318 条。余数或距离相同时优先输入位置靠前的类别或记录。

安装模型依赖：`python -m pip install "dnakit[neural]"`。
以下示例假设 `candidates.csv` 包含 `sequence`、`label` 列，其他列会随选中记录保留：

```python
from dnakit.datasets import DiversitySamplingConfig, diversity_sample
from dnakit.representations import RepresentationConfig

result = diversity_sample(
    "candidates.csv",
    n_samples=400,
    diversity="high",
    output_path="diversity_400",
    config=DiversitySamplingConfig(
        label_key="label",
        seed=42,
        representation=RepresentationConfig(
            model="dnabert2",
            pooling="mean",
            allow_remote_code=True,  # 使用官方 checkpoint 中的模型代码
            # checkpoint_path="/path/to/dnabert2-117m",  # 可复用本地模型
        ),
    ),
)
print(len(result.records))
for allocation in result.allocations:
    print(allocation.label, allocation.input_count, allocation.n_samples)
```

低多样性模式的调用方式：

```python
low_result = diversity_sample(
    "candidates.csv",
    n_samples=400,
    diversity="low",
    low_n_clusters=5,
    output_path="low_diversity_400",
    config=DiversitySamplingConfig(
        seed=42,
        representation=RepresentationConfig(model="dnabert2", allow_remote_code=True),
    ),
)
```

低多样性模式先在每个标签内部聚类。对于每个簇中心，从**该标签全部候选**中取配额数量的
最近邻，得到最多 5 个候选组。计算每组记录到**该组平均 embedding** 的平均欧氏距离，
选距离最小的一组；评分相同时按排序后的原始行索引打破平局。候选组可跨越簇边界，
因此单个簇不足额时也能选满。这个过程近似寻找低多样性子集，不保证全局最优。

文件输入会写出 `diversity_400/diversity_sample.csv` 和 `diversity_400/result.json`。
CSV 保留原表列及所选记录的输入顺序；JSON 包含名额分配、原始行索引、簇成员、中心及模型信息。
也可直接传入 `Iterable[DNARecord]`，将类别放在 `record.metadata["label"]`，从
`result.records` 获取保留序列、标签、注释和元数据的原始记录。

- `selected_indices` 是相对于本次输入的零起始索引；`cluster_labels` 是每条输入对应的簇编号。
  低多样性模式下，最终选择以 `records` 和 `selected_indices` 为准，不能用每簇代表拼出采样结果。
  结果同时记录 `diversity`、`low_n_clusters` 和 `sampling_strategy`；实际簇数见 `len(result.clusters)`。
- `label_key=None` 表示对整个无标签候选池采样。设置标签字段时，每条记录必须有非 null 的
  字符串、数字或布尔标签；不同类型的标签分别处理，如 `1` 与 `"1"`。
- 比例分配可能使极小类别获得 0 个名额；该类记录的 `cluster_labels` 为 `None`。
- 默认 `normalize=True`；设为 `False` 可直接在原始 embedding 的欧氏空间聚类。
  默认不做 PCA。其他 DNAKit 模型可通过 `RepresentationConfig(model=...)` 选择。
- `n_samples` 必须为正整数，不能超过候选数。`high` 模式特征不足以形成指定数量的簇时
  明确报错；`low` 模式将每个标签的簇数限制在其不同 embedding 数以内，并允许选取
  embedding 相同的不同输入行。两种模式均不重复采样同一行。拟合出现空簇时明确报错。
  默认输入上限为 10,000 条，同时受采样配置和表征配置的 `max_records` 限制。
- 本方法在模型特征空间中选择代表性或紧密子集，不保证最大或最小两两序列距离，也不做跨标签序列去重。
  同一 seed 的复现还取决于输入顺序、模型权重、依赖版本和计算环境。

验证覆盖两种模式的合成特征 840 → 400 采样、82/318 配额、紧密区域选择、重复特征和 CSV 导出；
不代表已对用户的候选序列完成 DNABERT2 推理或验证生物学多样性提升。
模型及算法参考：[DNABERT-2 官方实现](https://github.com/MAGICS-LAB/DNABERT_2)、
[scikit-learn KMeans](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html)。
