# 去重

对DNA序列进行去重，从而得到去重之后的DNA序列。

## 1) 标准去重 {#1-data-001}

- **作用：** 按序列逐字符完全相等进行分组，每组只保留一条代表记录，同时返回被合并记录及其代表映射，用于无损删除精确副本。
- **计算方法：** 以完全相同的序列作为分组键，按配置选择代表并处理 metadata 冲突。

- **API：** `dnakit.datasets.deduplicate(records[必须], equivalence[可选], config[可选])`；`config` 使用 `dnakit.datasets.DeduplicationConfig`。

- **输入：** 必填 `DNARecord` 集合；可选代表策略、冲突字段/策略和 metadata 合并。
- **输出：** 返回 `DeduplicationResult`，包含去重后记录 `records`、等价组 `groups` 和去重数量。

- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate

records = [
    DNARecord(DNASequence("ACGT"), "a"),
    DNARecord(DNASequence("TTAA"), "b"),
    DNARecord(DNASequence("ACGT"), "c"),
]
result = deduplicate(records, equivalence="exact")
print(result.records.ids, result.groups[0].member_ids)
```

- **示例结果：**

```text
('a', 'b') ('a', 'c')
```

## 2) 反向互补去重 {#2-data-002}

- **作用：** 比较序列的正向与反向互补形式，把仅链方向不同但内容相同的记录归为一组，避免双链方向造成重复计数。
- **计算方法：** 将序列与反向互补形式归入同一等价组，再选择代表。

- **API：** `dnakit.datasets.deduplicate(records[必须], equivalence[可选], config[可选])`；本项设置 `equivalence="reverse_complement"`。

- **输入：** 必填 `DNARecord` 集合；可选代表选择和 metadata 冲突配置。
- **输出：** 返回 `DeduplicationResult`，包含去重后记录 `records`、等价组 `groups` 和去重数量。

- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate

records = [
    DNARecord(DNASequence("AAGC"), "forward"),
    DNARecord(DNASequence("GCTT"), "reverse"),
]
result = deduplicate(records, equivalence="reverse_complement")
print(result.groups[0].member_ids, result.groups[0].orientations)
```

- **示例结果：**

```text
('forward', 'reverse') ('forward', 'reverse_complement')
```

## 3) 环状等价去重 {#3-data-003}

- **作用：** 对环状序列进行旋转等价比较，把存储起点不同但环上碱基顺序相同的记录合并，并保留代表记录和成员关系。
- **计算方法：** 将循环旋转规范化为等价键；可选同时将反向互补纳入等价关系。

- **API：** `dnakit.datasets.deduplicate(records[必须], equivalence[可选], config[可选])`；本项使用 `equivalence="circular"` 或 `equivalence="circular_reverse_complement"`。

- **输入：** 必填显式标记 `topology="circular"` 的记录；可选是否同时考虑反向互补。
- **输出：** 返回 `DeduplicationResult`，包含去重后记录 `records`、等价组 `groups` 和去重数量。

- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate

records = [
    DNARecord(DNASequence("AACG", topology="circular"), "origin-0"),
    DNARecord(DNASequence("CGAA", topology="circular"), "origin-2"),
]
result = deduplicate(records, equivalence="circular")
print(result.records.ids, result.groups[0].rotation_offsets)
```

- **示例结果：**

```text
('origin-0',) (0, 2)
```

## 4) IUPAC-aware去重 {#4-data-004-iupac-aware}

- **作用：** 把 IUPAC 模糊字符解释为碱基集合，逐位判断两条序列是完全相同、存在兼容可能还是明确冲突，用于含模糊碱基数据的等价检查。
- **计算方法：** 逐位比较 IUPAC 所代表的碱基集合，并按相应等价和兼容规则分组。

- **API：** `dnakit.datasets.deduplicate_iupac(records[必须], config[可选])`；`config` 使用 `dnakit.datasets.IUPACDeduplicationConfig`。

- **输入：** 必填可含 IUPAC 字符的记录；可选代表策略和两两比较上限。
- **输出：** 返回 `IUPACDeduplicationResult`，包含等价组及组内关系。

- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import deduplicate_iupac

records = [
    DNARecord(DNASequence("A", alphabet="iupac"), "a"),
    DNARecord(DNASequence("N", alphabet="iupac"), "n"),
    DNARecord(DNASequence("G", alphabet="iupac"), "g"),
]
result = deduplicate_iupac(records)
print([(group.member_ids, group.relation) for group in result.groups])
```

- **示例结果：**

```text
[(('a', 'n'), 'compatible'), (('g',), 'singleton')]
```

## 5) 近似去重 {#5-data-005}

- **作用：** 使用 Identity、Edit distance 或 k-mer 相似度比较序列，把超过阈值的近似副本合并成组，并返回代表记录和分组依据。
- **计算方法：** 计算成对相似度并按阈值形成簇，再为每簇保留代表记录。

- **API：** `dnakit.datasets.deduplicate_approximate(records[必须], config[必须])`；`config` 使用 `dnakit.datasets.ClusterConfig`。

- **输入：** 必填记录集合和相似度方法/阈值；可选 k、canonical、代表策略和资源上限。
- **输出：** 返回 `ClusteringResult`，包含近似重复分组和代表记录。

- **示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.datasets import ClusterConfig, deduplicate_approximate

records = [
    DNARecord(DNASequence("AAAA"), "a"),
    DNARecord(DNASequence("AAAT"), "b"),
    DNARecord(DNASequence("CCCC"), "c"),
]
result = deduplicate_approximate(
    records, config=ClusterConfig(method="identity", threshold=0.7)
)
print(result.labels, result.representatives.ids)
```

- **示例结果：**

```text
(0, 0, 1) ('a', 'c')
```
