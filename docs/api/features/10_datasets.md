# 数据划分

按照随机、分层、相似度或元数据约束划分 DNA 数据集，从而得到训练、验证和测试等相互隔离的数据子集。

所有随机或启发式过程都记录 seed、配置、分组和资源上限。

## 1) DATA-012 · 随机与稳定哈希划分

**作用：** 按目标比例把记录随机分配到 train、valid、test，并结合稳定 ID 与 seed 生成可复现结果，返回各子集和分配清单。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="random"` 或 `method="hash"`。

**输入：** 必填记录和总和为 1 的比例；可选 seed、shuffle 和是否保持子集内原顺序。`hash` 模式要求每条记录具有稳定且唯一的 `record.id`。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
result = split(
    records,
    config=SplitConfig(
        method="random", ratios={"train": 0.8, "test": 0.2}, seed=17
    ),
)
print(dict(result.counts))
```

**示例结果：**

```text
{'train': 8, 'test': 2}
```

<span id="1-hash"></span>**顺序无关的 `hash` 划分**

`hash` 模式在 `shuffle=True` 时使用版本化 SHA-256 计算 `seed + record.id` 的稳定排序键，在 `shuffle=False` 时按 `record.id` 的稳定字节序排列，再按比例分配记录。因此，同一批记录即使输入顺序变化，按 `record.id` 查看时分类仍完全一致；它不使用 Python 内置的进程随机化 `hash()`。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import SplitConfig, split

records = [DNARecord(DNASequence("AC"), f"r{i}") for i in range(10)]
config = SplitConfig(
    method="hash",
    ratios={"train": 0.6, "valid": 0.2, "test": 0.2},
    seed=17,
    preserve_order=False,
)
first = split(records, config=config)
second = split(list(reversed(records)), config=config)
first_by_id = {item.record_id: item.split for item in first.assignments}
second_by_id = {item.record_id: item.split for item in second.assignments}
print(first_by_id == second_by_id)
print(first.get("train").ids == second.get("train").ids)
```

**示例结果：**

```text
True
True
```

## 2) DATA-013 · 分层随机划分

**作用：** 在 train、valid、test 中尽量保持指定类别标签的原始比例，返回分层后的子集及标签统计，减少类别分布偏移。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="stratified"`。

**输入：** 必填记录、比例和 `metadata_key`；可选 seed、缺失字段策略和顺序策略。

**示例代码：**

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

**示例结果：**

```text
{'train': 4, 'test': 4}
```

## 3) DATA-014 · 相似度划分

**作用：** 先按相似度阈值建立序列组，再以整组为单位分配 train、valid、test，防止近似序列跨集合形成数据泄漏。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="similarity"`。

**输入：** 必填记录和比例；可选 k、阈值、IUPAC/Gap 策略、seed 和规模上限。

**示例代码：**

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

**示例结果：**

```text
{'a': 'train', 'b': 'train', 'c': 'test', 'd': 'test'}
```

## 4) DATA-015 · Cluster split

**作用：** 使用已有 cluster 标签作为不可拆分分组进行数据划分，保证同一 cluster 的全部记录只进入一个子集，并报告实际比例偏差。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.SplitConfig`，本项设置 `method="group", metadata_key="cluster"`。

**输入：** 必填记录、每条记录的 cluster metadata 和比例；可选 seed、缺失值策略。

**示例代码：**

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

**示例结果：**

```text
{'a1': 'train', 'a2': 'train', 'b1': 'test'}
```

## 5) DATA-016 · 物种划分

**作用：** 按物种 metadata 将记录成组后整体划分，确保同一物种不同时出现在训练集和评估集，用于检验跨物种泛化。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。

**输入：** 必填 species metadata key 和比例；可选 seed、缺失值策略。

**示例代码：**

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

**示例结果：**

```text
{'h1': 'train', 'h2': 'train', 'm1': 'test'}
```

## 6) DATA-017 · 染色体划分

**作用：** 按染色体 metadata 将记录成组后整体划分，避免同一染色体的高度相关区域跨 train、valid、test 造成泄漏。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。

**输入：** 必填 chromosome metadata key 和比例；可选 seed、缺失值策略。

**示例代码：**

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

**示例结果：**

```text
{'c1-a': 'train', 'c1-b': 'train', 'c2-a': 'test'}
```

## 7) DATA-018 · 个体划分

**作用：** 按个体或 donor metadata 对样本分组并整体划分，确保同一个体不跨数据子集，用于避免个体特异信息泄漏。


**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig` 的 `group` 模式。

**输入：** 必填 individual metadata key 和比例；可选 seed、缺失值策略。

**示例代码：**

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

**示例结果：**

```text
{'p1-a': 'train', 'p1-b': 'train', 'p2-a': 'test'}
```

## 8) 按自定义 label 划分

**作用：** 按调用方指定的任意 metadata label（如 family、locus、batch）对记录分组，确保相同 label 值只进入一个子集，并返回分组分配结果。

**API：** `dnakit.datasets.split(records[必须], config[必须])`；使用 `dnakit.datasets.SplitConfig(method="group", metadata_key=label)`。

**输入：** 必填记录、划分比例和自定义 `label`。`label` 是每条记录中的 metadata 字段名，例如 `donor`、`family`、`locus`、`batch` 或其他调用方提供的字段。

**示例代码：**

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

**示例结果：**

```text
{'f1-a': 'train', 'f1-b': 'train', 'f2-a': 'test'}
```

## 9) DATA-023 · 泄漏检测

**作用：** 对 train、valid、test 等集合执行跨集合精确或近似比较，返回泄漏序列对、相似度、涉及子集和汇总比例，用于验证划分独立性。


**API：** `dnakit.datasets.detect_leakage(splits[必须], config[可选])`；`config` 使用 `dnakit.datasets.LeakageConfig`。

**输入：** 必填 split 名称到 `DNASet` 的映射；可选 identity/edit/k-mer/fingerprint、阈值和资源上限。

**示例代码：**

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

**示例结果：**

```text
True 1
```
