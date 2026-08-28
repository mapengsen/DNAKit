# 通用搜索

在单条或多条 DNA 序列中进行精确、子序列、近似、反向互补和最近邻搜索。

## 1) SIM-001 · Exact search

**作用：** 在记录集合中查找与查询序列逐字符完全相同的记录，返回目标索引、起止坐标和链方向，用于成员检查、精确副本定位和小规模参考查询。

**API：** `dnakit.similarity.exact_search(query[必须], targets[必须], reverse_complement[可选], merge_strands[可选], max_targets[可选], max_matches[可选])`。

**输入：** 必填查询序列和目标序列/`DNASet`；可选反向互补搜索、链合并及资源上限。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import exact_search

result = exact_search(
    DNASequence("AC"),
    [DNASequence("AC"), DNASequence("TACG")],
)
print([(hit.target_index, hit.start, hit.end) for hit in result.matches])
```

**示例结果：**

```text
[(0, 0, 2)]
```


## 2) SIM-002 · Subsequence search

**作用：** 在目标长序列中查找短查询序列的全部精确出现位置，返回起止坐标和匹配方向，用于定位 motif、引物或已知片段。

**API：** `dnakit.similarity.subsequence_search(query[必须], target[必须], strand[可选], overlapping[可选], merge_strands[可选], max_targets[可选], max_matches[可选])`。

**输入：** 必填 query 和 target；可选 `strand`、`overlapping`、链合并及资源上限。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import subsequence_search

result = subsequence_search(DNASequence("ANA", alphabet="iupac"),
                            DNASequence("ANANA", alphabet="iupac"))
print([(hit.start, hit.end) for hit in result.matches])
```

**示例结果：**

```text
[(0, 3), (2, 5)]
```


## 3) SIM-004 · Approximate matching

**作用：** 在目标序列中查找允许指定数量替换、插入或删除的近似命中，返回坐标、距离和匹配片段，用于容忍突变或测序误差的搜索。

**API：** `dnakit.similarity.approximate_search(query[必须], targets[必须], max_distance[必须], substitution_cost[可选], insertion_cost[可选], deletion_cost[可选], reverse_complement[可选], max_targets[可选], max_matches[可选], max_cells[可选])`。

**输入：** 必填 query、target 和 `max_distance`；可选三类编辑代价、反向互补和资源上限。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import approximate_search

result = approximate_search(DNASequence("ACG"), DNASequence("TTACGATG"),
                            max_distance=1)
print(any(hit.start == 2 and hit.distance == 0 for hit in result.matches))
```

**示例结果：**

```text
True
```

## 4) SIM-005 · 反向互补搜索

**作用：** 同时使用查询序列及其反向互补序列进行搜索，返回每个命中的坐标和链方向，避免遗漏位于反向链的相同片段。

**API：** `dnakit.similarity.reverse_complement_search(query[必须], target[必须], overlapping[可选], merge_strands[可选], max_targets[可选], max_matches[可选])`。

**输入：** 必填 query 和 target；可选是否允许重叠、是否合并回文重复命中。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import reverse_complement_search

result = reverse_complement_search(DNASequence("ATG"), DNASequence("GGCATCC"))
print([(hit.start, hit.end, hit.strand.value) for hit in result.matches])
```

**示例结果：**

```text
[(2, 5, 'reverse')]
```

## 5) SIM-014 · 最近邻搜索

**作用：** 使用查询序列的 k-mer Sketch 在已有索引中筛选相似度最高的 Top-k 记录，返回记录 ID、排名和近似相似度，适合快速候选检索。

**API：** `dnakit.similarity.build_sketch_index(records[必须], k[可选], num_hashes[可选], canonical[可选], seed[可选], max_records[可选])`、`dnakit.similarity.nearest_neighbors(query[必须], index[必须], top_k[可选], min_similarity[可选])`。

**输入：** 必填 query 和 `SketchIndex`；可选 `top_k`、最低相似度及构建索引时的 sketch 参数。

**示例代码：**

```python
from dnakit import DNARecord, DNASequence
from dnakit.similarity import build_sketch_index, nearest_neighbors

records = [
    DNARecord(DNASequence("ACGTACGT"), "a"),
    DNARecord(DNASequence("ACGTTCGT"), "b"),
]
index = build_sketch_index(records, k=3, num_hashes=100)
result = nearest_neighbors(records[0], index, top_k=2)
print([hit.record_id for hit in result.hits])
```

**示例结果：**

```text
['a', 'b']
```


## 6) SIM-015 · 数据库索引

**作用：** 为参考序列集合创建带参数和记录 ID 的可复用 Sketch 索引，并支持保存、校验和重新加载，避免每次 Top-k 查询都重复计算参考摘要。

**API：** `dnakit.similarity.build_sketch_index(records[必须], k[可选], num_hashes[可选], canonical[可选], seed[可选], max_records[可选])`、`dnakit.similarity.save_sketch_index(index[必须], path[必须], overwrite[可选])`、`dnakit.similarity.load_sketch_index(path[必须])`。

**输入：** 必填 `DNARecord` 集合；可选 k、hash 数、canonical、seed、保存路径和记录上限。

**示例代码：**

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from dnakit import DNARecord, DNASequence
from dnakit.similarity import build_sketch_index, load_sketch_index, save_sketch_index

records = [DNARecord(DNASequence("ACGT"), "r1")]
with TemporaryDirectory() as directory:
    path = Path(directory) / "sketch-index.json"
    index = build_sketch_index(records, k=2, num_hashes=20)
    digest = save_sketch_index(index, path)
    print(load_sketch_index(path) == index, len(digest))
```

**示例结果：**

```text
True 64
```
