# 数据划分

按照随机、分层、相似度或元数据约束划分 DNA 数据集，从而得到训练、验证和测试等相互隔离的数据子集。

所有随机或启发式过程都记录 seed、配置、分组和资源上限。

这些记录集合接口可直接接收 `DNA(["ACGT", "TTAA", ...])`：`DNA` 会按顺序提供 `DNARecord`。需要稳定 ID 或每条记录的 metadata 时，改用 `DNA([{"id": "r1", "sequence": "ACGT", "metadata": {...}}, ...])`，也可继续传 `DNARecord` 迭代器。

## 1) 随机与稳定哈希划分 {#1-data-012}

- **作用：** 按目标比例把记录分配到 `train / valid / test` 等命名子集，或用固定 seed 生成 K-fold 划分并选择一折作为 `valid`、其余折作为 `train`。返回各子集和分配清单。
- **处理方法：** 随机模式按 seed 打乱后按比例分配；hash 模式按稳定 ID 排序；K-fold 固定划分后选择验证折。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项支持 `method="random"`、`method="hash"` 和 `method="kfold"`。
- **输入：** 必填记录；比例划分通过 `ratios` 指定总和为 1 的比例，默认 `train=0.8, valid=0.1, test=0.1`。K-fold 需要 scikit-learn，通过 `n_splits` 指定折数、`fold_id` 指定验证折。默认 `seed=42`、`shuffle=True`。`hash` 模式要求每条记录具有稳定且唯一的 `record.id`。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.datasets import SplitConfig, split

records = DNA([{"id": f"r{i}", "sequence": "AC"} for i in range(10)])
result = split(
    records,
    config=SplitConfig(
        method="random", ratios={"train": 0.8, "test": 0.2}, seed=17
    ),
)
print(dict(result.counts))
```

- **示例结果：**

```text
{'train': 8, 'test': 2}
```

<span id="1-hash"></span>稳定哈希、K-fold 示例、依赖安装与旧版 seed 兼容说明见[FAQ：详细说明](../../faq/tutorial_integrations.md#dataset-splits)。

正反链共享 `pair_id` 时，可把它放在每条记录的 metadata 中，使用 `SplitConfig(method="group", metadata_key="pair_id")` 按组划分。文件形式的数据可直接使用 `export_grouped_split`，它对每个 `pair_id` 只分配一次，自动写出分配表和各子集 CSV，无需为每一行构造 `DNARecord`。`method="kfold"` 按输入记录逐条分折；若要让配对分块同折，应先为每个 `pair_id` 建立一条代表记录，再将折号映射回两条链。Fig.4 的 A05 对 hg38 分块使用前述文件接口，而 `--n-splits` 控制两个基准数据集官方训练集的折数。

### 文件型分组划分与 K-fold 导出

`dnakit.datasets.export_grouped_split(source, output_dir, group_field=..., config=...)` 按组字段进行随机比例划分。默认 `SplitConfig()` 是种子 42 的 8:1:1 划分。它在输出目录生成 `assignment.csv` 和 `train.csv`、`valid.csv`、`test.csv`；各组内所有行进入同一子集，行顺序保持原始顺序。`sample_sizes` 可按子集指定蓄水池抽样数量，`sample_filter` 控制候选行，`sample_digest_field` 为抽样结果添加 `<字段>_sha256` 列；抽样 CSV 写在 `analysis/`。返回值包含 `seed`、`group_counts` 和 `analysis_counts`。当前文件接口支持 `SplitConfig(method="random")`，输入需提供可读取的行文件，如 CSV 或 Parquet；输出为 CSV。

```python
from dnakit.datasets import SplitConfig, export_grouped_split, export_kfold_table

summary = export_grouped_split(
    "chunks.csv", "split_output", group_field="pair_id",
    config=SplitConfig(seed=42),
    sample_sizes={"train": 800, "valid": 100, "test": 100},
    sample_filter=lambda row: row["strand"] == "+",
    sample_digest_field="sequence",
    overwrite=True,
)
fold_counts = export_kfold_table(
    "benchmark_train.csv", "benchmark_train_with_folds.csv",
    config=SplitConfig(method="kfold", n_splits=10, seed=42),
    label_field="label", overwrite=True,
)
```

`export_kfold_table` 保留源文件行顺序，附加从 0 开始的 `fold` 列；传入 `label_field` 时返回每折训练和验证标签计数，可用 `dnakit.io.save_json` 保存。默认 10 折、种子 42，需要 scikit-learn。它导出所有折的编号，因此 `SplitConfig.fold_id` 不用于选择单独一折。

## 2) 分层随机划分 {#2-data-013}

- **作用：** 在 train、valid、test 中尽量保持指定类别标签的原始比例，返回分层后的子集及标签统计，减少类别分布偏移。
- **处理方法：** 按标签分组，再按目标比例分配各组记录。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="stratified"`。
- **输入：** 必填记录、比例和 `metadata_key`；可选 seed、缺失字段策略和顺序策略。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), f"p{i}", metadata={"label": "p"})
    for i in range(4)
] + [
    DNARecord(DNASequence("GT"), f"n{i}", metadata={"label": "n"})
    for i in range(4)
]
result = split(
    records,
    config=SplitConfig(
        method="stratified", ratios={"train": 0.5, "test": 0.5},
        metadata_key="label", seed=2,
    ),
)
print(dict(result.counts))
```

- **示例结果：**

```text
{'train': 4, 'test': 4}
```

## 3) 相似度划分 {#3-data-014}

- **作用：** 先按相似度阈值建立序列组，再以整组为单位分配 train、valid、test，防止近似序列跨集合形成数据泄漏。
- **处理方法：** 按 k-mer 相似度阈值建立组，再将整组分配到子集。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="similarity"`。
- **输入：** 必填记录和比例；可选 k、阈值、IUPAC/Gap 策略、seed 和规模上限。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
    DNARecord(DNASequence("CCCG"), "d"),
]
result = split(
    records,
    config=SplitConfig(
        method="similarity", ratios={"train": 0.5, "test": 0.5},
        similarity_k=2, similarity_threshold=0.5, seed=5,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

- **示例结果：**

```text
{'a': 'train', 'b': 'train', 'c': 'test', 'd': 'test'}
```

## 4) Cluster split {#4-data-015-cluster-split}

- **作用：** 使用已有 cluster 标签作为不可拆分分组进行数据划分，保证同一 cluster 的全部记录只进入一个子集，并报告实际比例偏差。
- **处理方法：** 按聚类字段形成组，再将每个完整簇分配到一个子集。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="group", metadata_key="cluster"`。
- **输入：** 必填记录、每条记录的 cluster metadata 和比例；可选 seed、缺失值策略。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "a1", metadata={"cluster": "A"}),
    DNARecord(DNASequence("AG"), "a2", metadata={"cluster": "A"}),
    DNARecord(DNASequence("GT"), "b1", metadata={"cluster": "B"}),
]
result = split(
    records,
    config=SplitConfig(
        method="group", ratios={"train": 0.5, "test": 0.5},
        metadata_key="cluster", seed=11,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

- **示例结果：**

```text
{'a1': 'train', 'a2': 'train', 'b1': 'test'}
```

## 5) 物种划分 {#5-data-016}

- **作用：** 按物种 metadata 将记录成组后整体划分，确保同一物种不同时出现在训练集和评估集，用于检验跨物种泛化。
- **处理方法：** 按物种元数据分组，同一物种只进入一个子集。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。
- **输入：** 必填 species metadata key 和比例；可选 seed、缺失值策略。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "h1", metadata={"species": "human"}),
    DNARecord(DNASequence("AG"), "h2", metadata={"species": "human"}),
    DNARecord(DNASequence("GT"), "m1", metadata={"species": "mouse"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="species"
))
print({item.record_id: item.split for item in result.assignments})
```

- **示例结果：**

```text
{'h1': 'train', 'h2': 'train', 'm1': 'test'}
```

## 6) 染色体划分 {#6-data-017}

- **作用：** 按染色体 metadata 将记录成组后整体划分，避免同一染色体的高度相关区域跨 train、valid、test 造成泄漏。
- **处理方法：** 按染色体元数据分组，同一染色体只进入一个子集。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。
- **输入：** 必填 chromosome metadata key 和比例；可选 seed、缺失值策略。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "c1-a", metadata={"chromosome": "chr1"}),
    DNARecord(DNASequence("AG"), "c1-b", metadata={"chromosome": "chr1"}),
    DNARecord(DNASequence("GT"), "c2-a", metadata={"chromosome": "chr2"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="chromosome"
))
print({item.record_id: item.split for item in result.assignments})
```

- **示例结果：**

```text
{'c1-a': 'train', 'c1-b': 'train', 'c2-a': 'test'}
```

## 7) 个体划分 {#7-data-018}

- **作用：** 按个体或 donor metadata 对样本分组并整体划分，确保同一个体不跨数据子集，用于避免个体特异信息泄漏。
- **处理方法：** 按个体元数据分组，同一个体只进入一个子集。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。
- **输入：** 必填 individual metadata key 和比例；可选 seed、缺失值策略。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [
    DNARecord(DNASequence("AC"), "p1-a", metadata={"individual": "p1"}),
    DNARecord(DNASequence("AG"), "p1-b", metadata={"individual": "p1"}),
    DNARecord(DNASequence("GT"), "p2-a", metadata={"individual": "p2"}),
]
result = split(records, config=SplitConfig(
    method="group", ratios={"train": 0.5, "test": 0.5}, metadata_key="individual"
))
print({item.record_id: item.split for item in result.assignments})
```

- **示例结果：**

```text
{'p1-a': 'train', 'p1-b': 'train', 'p2-a': 'test'}
```

## 8) 按自定义 label 划分

- **作用：** 按调用方指定的任意 metadata label（如 family、locus、batch）对记录分组，确保相同 label 值只进入一个子集，并返回分组分配结果。
- **处理方法：** 按指定 metadata 字段的 label 值分组，再将整组分配到子集。
- **API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig(method="group", metadata_key=label)`。
- **输入：** 必填记录、划分比例和自定义 `label`。`label` 是每条记录中的 metadata 字段名，例如 `donor`、`family`、`locus`、`batch` 或其他调用方提供的字段。
- **输出：** 返回 `SplitResult`，包含命名子集 `subsets`、各子集数量 `counts` 和逐条分配 `assignments`。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

label = "family"  # 也可以替换为 donor、locus 或其他 metadata 字段名
records = [
    DNARecord(DNASequence("AC"), "f1-a", metadata={label: "f1"}),
    DNARecord(DNASequence("AG"), "f1-b", metadata={label: "f1"}),
    DNARecord(DNASequence("GT"), "f2-a", metadata={label: "f2"}),
]
result = split(
    records,
    config=SplitConfig(
        method="group",
        ratios={"train": 0.5, "test": 0.5},
        metadata_key=label,
    ),
)
print({item.record_id: item.split for item in result.assignments})
```

- **示例结果：**

```text
{'f1-a': 'train', 'f1-b': 'train', 'f2-a': 'test'}
```

## 9) 泄漏检测 {#9-data-023}

- **作用：** 对 train、valid、test 等集合执行跨集合精确或近似比较，返回泄漏序列对、相似度、涉及子集和汇总比例，用于验证划分独立性。
- **处理方法：** 比较不同子集的序列或分组信息，报告跨集合重复或相似命中。
- **API：** `dnakit.datasets.detect_leakage(splits[必须], config[可选])`；`config` 使用 `dnakit.datasets.LeakageConfig`。
- **输入：** 必填 split 名称到 `DNASet` 的映射；可选 identity/edit/k-mer/fingerprint、阈值和资源上限。
- **输出：** 返回 `LeakageReport`，包含是否泄漏、泄漏事件数和跨子集命中。
- **示例代码：**

```python
from dnakit import DNARecord, DNASequence, DNASet
from dnakit.datasets import LeakageConfig, detect_leakage

splits = {
    "train": DNASet([DNARecord(DNASequence("AAAA"), "a")]),
    "test": DNASet([DNARecord(DNASequence("AAAA"), "same")]),
}
report = detect_leakage(
    splits, config=LeakageConfig(method="identity", threshold=0.9)
)
print(report.has_leakage, report.exact_event_count)
```

- **示例结果：**

```text
True 1
```

<!-- BEGIN GENERATED API: encode-collection -->

<span id="encode-collection"></span>

## 10) 比较 ENCODE 实验的联合分析兼容性 {#dnakit.encode.compare_experiments}

- **作用：** 比较两个已追踪实验的兼容性。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.compare_experiments(accession1, accession2, *, verbose=True)`
- **输入：** 必须提供 `accession1`, `accession2`。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，包含兼容性结论、问题、警告和建议。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.compare_experiments(
    accession1='ENCSR133RZO',
    accession2='ENCSR000AKS',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含兼容性结论、问题、警告和建议。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.compare_experiments-parameters)

## 11) 将追踪的 ENCODE 实验导出为表格 {#dnakit.encode.export_data}

- **作用：** 把本地追踪的实验导出为 CSV、TSV 或 JSON 文本。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.export_data(format='csv', assay_title=None, organism=None, organ=None, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回指定 CSV、TSV 或 JSON 格式的表格文本。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.export_data(
    format='csv',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回指定 CSV、TSV 或 JSON 格式的表格文本。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.export_data-parameters)

## 12) 获取追踪实验的论文与参考文献 {#dnakit.encode.get_citations}

- **作用：** 获取追踪实验的论文引用，可导出 BibTeX 或 RIS。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.get_citations(accession=None, export_format='json', *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回指定格式的文献引用文本，例如 BibTeX 或 RIS。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.get_citations(
    export_format='bibtex',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回指定格式的文献引用文本，例如 BibTeX 或 RIS。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.get_citations-parameters)

## 13) 查询衍生文件的来源记录 {#dnakit.encode.get_provenance}

- **作用：** 查询衍生文件的来源链。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.get_provenance(file_path=None, source_accession=None, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，包含来源链或衍生文件列表。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.get_provenance(
    file_path='analysis/consensus_peaks.bed',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含来源链或衍生文件列表。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.get_provenance-parameters)

## 14) 查询实验关联的外部参考标识 {#dnakit.encode.get_references}

- **作用：** 查询实验关联的外部引用。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.get_references(experiment_accession=None, reference_type=None, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，列出关联的外部引用。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.get_references(
    reference_type='doi',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，列出关联的外部引用。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.get_references-parameters)

## 15) 为追踪实验关联外部参考标识 {#dnakit.encode.link_reference}

- **作用：** 关联实验与 DOI、PMID、GEO 等外部标识。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.link_reference(experiment_accession, reference_type, reference_id, description='', *, verbose=True)`
- **输入：** 必须提供 `experiment_accession`, `reference_type`, `reference_id`。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，报告引用关联结果。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.link_reference(
    experiment_accession='ENCSR133RZO',
    reference_type='doi',
    reference_id='10.1038/nature11247',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，报告引用关联结果。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.link_reference-parameters)

## 16) 列出本地追踪的 ENCODE 实验 {#dnakit.encode.list_tracked}

- **作用：** 列出本地追踪的实验。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.list_tracked(assay_title=None, organism=None, organ=None, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，包含已追踪实验的元数据表和追踪统计。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.list_tracked(
    assay_title='ATAC-seq',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含已追踪实验的元数据表和追踪统计。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.list_tracked-parameters)

## 17) 登记衍生文件与数据来源 {#dnakit.encode.log_derived_file}

- **作用：** 登记衍生文件及其来源、工具和参数。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.log_derived_file(file_path, source_accessions, description='', file_type='', tool_used='', parameters='', *, verbose=True)`
- **输入：** 必须提供 `file_path`, `source_accessions`。 先用文件查询功能取得实际编号，替换 `ENCFF_REPLACE_ME`。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，包含新建来源记录的 ID。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.log_derived_file(
    file_path='analysis/consensus_peaks.bed',
    source_accessions=['ENCFF_REPLACE_ME'],
    description='Consensus peak set',
    file_type='bed',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含新建来源记录的 ID。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.log_derived_file-parameters)

## 18) 按字段汇总追踪的实验集合 {#dnakit.encode.summarize_collection}

- **作用：** 按实验、组织等字段汇总本地实验集合。
- **处理方法：** 读取或更新本地 ENCODE 追踪记录，按所选操作整理结果。
- **API：** `dnakit.encode.summarize_collection(assay_title=None, organism=None, organ=None, *, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，按多个维度汇总实验数量。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.summarize_collection(
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，按多个维度汇总实验数量。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.summarize_collection-parameters)

## 19) 在本地追踪 ENCODE 实验及其来源信息 {#dnakit.encode.track_experiment}

- **作用：** 把实验加入本地追踪目录。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.track_experiment(accession, fetch_publications=True, fetch_pipelines=True, notes='', *, verbose=True)`
- **输入：** 必须提供 `accession`。 集合操作使用本地追踪记录；比较、关联引用和导出前，先用 `track_experiment` 登记实验。文件来源查询还需先登记衍生文件。
- **输出：** 返回 JSON 文本，包含追踪结果、找到的文献与分析流程信息。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.track_experiment(
    accession='ENCSR133RZO',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含追踪结果、找到的文献与分析流程信息。

[完整参数与使用说明](../../faq/api/10_datasets.md#dnakit.encode.track_experiment-parameters)

<!-- END GENERATED API: encode-collection -->
