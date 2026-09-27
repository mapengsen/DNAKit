# 通用搜索

在单条或多条 DNA 序列中进行精确、子序列、近似、反向互补和最近邻搜索。

## 1) Exact search {#1-sim-001-exact-search}

- **作用：** 在记录集合中查找与查询序列逐字符完全相同的记录，返回目标索引、起止坐标和链方向，用于成员检查、精确副本定位和小规模参考查询。
- **处理方法：** 将查询与每个目标的完整符号串逐字符比较。
- **API：** `dnakit.similarity.exact_search(query[必须], targets[必须], reverse_complement[可选], merge_strands[可选], max_targets[可选], max_matches[可选])`。
- **输入：** 必填查询序列和目标序列/`DNASet`；可选反向互补搜索、链合并及资源上限。
- **输出：** 返回 `SearchResult`，匹配位置与命中信息在 `matches`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import exact_search

result = exact_search(
    DNASequence("AC"),
    [DNASequence("AC"), DNASequence("TACG")],
)
print([(hit.target_index, hit.start, hit.end) for hit in result.matches])
```

- **示例结果：**

```text
[(0, 0, 2)]
```


## 2) Subsequence search {#2-sim-002-subsequence-search}

- **作用：** 在目标长序列中查找短查询序列的全部精确出现位置，返回起止坐标和匹配方向，用于定位 motif、引物或已知片段。
- **处理方法：** 扫描目标中的精确子串位置，按配置保留重叠命中。
- **API：** `dnakit.similarity.subsequence_search(query[必须], target[必须], strand[可选], overlapping[可选], merge_strands[可选], max_targets[可选], max_matches[可选])`。
- **输入：** 必填 query 和 target；可选 `strand`、`overlapping`、链合并及资源上限。
- **输出：** 返回 `SearchResult`，匹配位置与命中信息在 `matches`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import subsequence_search

result = subsequence_search(DNASequence("ANA", alphabet="iupac"),
                            DNASequence("ANANA", alphabet="iupac"))
print([(hit.start, hit.end) for hit in result.matches])
```

- **示例结果：**

```text
[(0, 3), (2, 5)]
```


## 3) Approximate matching {#3-sim-004-approximate-matching}

- **作用：** 在目标序列中查找允许指定数量替换、插入或删除的近似命中，返回坐标、距离和匹配片段，用于容忍突变或测序误差的搜索。
- **处理方法：** 用替换、插入和删除的加权编辑距离筛选阈值内命中。
- **API：** `dnakit.similarity.approximate_search(query[必须], targets[必须], max_distance[必须], substitution_cost[可选], insertion_cost[可选], deletion_cost[可选], reverse_complement[可选], max_targets[可选], max_matches[可选], max_cells[可选])`。
- **输入：** 必填 query、target 和 `max_distance`；可选三类编辑代价、反向互补和资源上限。
- **输出：** 返回 `ApproximateSearchResult`，示例展示 `matches`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import approximate_search

result = approximate_search(DNASequence("ACG"), DNASequence("TTACGATG"),
                            max_distance=1)
print(any(hit.start == 2 and hit.distance == 0 for hit in result.matches))
```

- **示例结果：**

```text
True
```

## 4) 反向互补搜索 {#4-sim-005}

- **作用：** 同时使用查询序列及其反向互补序列进行搜索，返回每个命中的坐标和链方向，避免遗漏位于反向链的相同片段。
- **处理方法：** 同时搜索查询及其反向互补序列，并记录链方向。
- **API：** `dnakit.similarity.reverse_complement_search(query[必须], target[必须], overlapping[可选], merge_strands[可选], max_targets[可选], max_matches[可选])`。
- **输入：** 必填 query 和 target；可选是否允许重叠、是否合并回文重复命中。
- **输出：** 返回 `SearchResult`，匹配位置与命中信息在 `matches`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import reverse_complement_search

result = reverse_complement_search(DNASequence("ATG"), DNASequence("GGCATCC"))
print([(hit.start, hit.end, hit.strand.value) for hit in result.matches])
```

- **示例结果：**

```text
[(2, 5, 'reverse')]
```

## 5) 最近邻搜索 {#5-sim-014}

- **作用：** 使用查询序列的 k-mer Sketch 在已有索引中筛选相似度最高的 Top-k 记录，返回记录 ID、排名和近似相似度，适合快速候选检索。
- **处理方法：** 比较查询与索引中的 k-mer Sketch，按近似相似度排序取 Top-k。
- **API：** `dnakit.similarity.build_sketch_index(records[必须], k[可选], num_hashes[可选], canonical[可选], seed[可选], max_records[可选])`、`dnakit.similarity.nearest_neighbors(query[必须], index[必须], top_k[可选], min_similarity[可选])`。
- **输入：** 必填 query 和 `SketchIndex`；可选 `top_k`、最低相似度及构建索引时的 sketch 参数。
- **输出：** 索引构建返回 `SketchIndex`；查询返回按相似度排序的最近邻命中。
- **示例代码：**

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

- **示例结果：**

```text
['a', 'b']
```


## 6) 数据库索引 {#6-sim-015}

- **作用：** 为参考序列集合创建带参数和记录 ID 的可复用 Sketch 索引，并支持保存、校验和重新加载，避免每次 Top-k 查询都重复计算参考摘要。
- **处理方法：** 为各记录计算带固定参数的 Sketch，将索引及完整性校验信息保存到 JSON。
- **API：** `dnakit.similarity.build_sketch_index(records[必须], k[可选], num_hashes[可选], canonical[可选], seed[可选], max_records[可选])`、`dnakit.similarity.save_sketch_index(index[必须], path[必须], overwrite[可选])`、`dnakit.similarity.load_sketch_index(path[必须])`。
- **输入：** 必填 `DNARecord` 集合；可选 k、hash 数、canonical、seed、保存路径和记录上限。
- **输出：** 返回 `SketchIndex`，保存索引记录及参数，可用于后续相似性查询。
- **示例代码：**

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

- **示例结果：**

```text
True 64
```

## 7) 建立本地 BLAST 数据库并解析命中 {#bio-4}

- **作用：** 在自己的核酸序列库中搜索查询序列。
- **计算方法：** 调用 `makeblastdb` 和 `blastn`，使用 SearchIO 解析 XML。
- **API：** `dnakit.bioinformatics.alignment.build_blast_database(fasta, output_directory)`; `dnakit.bioinformatics.alignment.local_blast(query, database, output_directory)`; `dnakit.bioinformatics.alignment.parse_search_results(path, format='blast-xml')`
- **输入：** 已有 `reference.fa`、`query.fa`；BLAST+ 在 PATH 中；需要 `bio-formats`。输出目录必须尚不存在。
- **输出：** 建库和搜索返回 `ToolResult`，结果 XML 路径在 `artifacts`；解析函数返回查询及其命中记录。
- **示例代码：**

```python
from dnakit.bioinformatics.alignment import build_blast_database, local_blast, parse_search_results

database = build_blast_database("reference.fa", "new_blast_db")
result = local_blast("query.fa", database.directory / "database", "new_blast_search")
for query in parse_search_results(result.artifacts["results.xml"]):
    print(query.id, len(query))
```

- **示例结果：**

返回 `ToolResult`，XML 路径在 `artifacts['results.xml']`；循环输出查询 ID 与命中数量。真实命中取决于输入，此处未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#alignment)。

## 8) 构建 Profile HMM 并搜索核酸序列 {#bio-5}

- **作用：** 从已有比对建立模型，再寻找相似序列片段。
- **计算方法：** 使用 `hmmbuild` 建模、`nhmmer` 搜索并解析核酸命中表；蛋白搜索可用 `hmmsearch`。
- **API：** `dnakit.bioinformatics.alignment.build_profile(alignment, output_directory, alphabet='dna')`; `dnakit.bioinformatics.alignment.search_profile(profile, sequences, output_directory, method='nhmmer')`; `dnakit.bioinformatics.alignment.parse_nhmmer_table(path)`
- **输入：** 已有 `aligned.sto` 和 `targets.fa`；HMMER 在 PATH 中；新输出目录。
- **输出：** 建模和搜索返回 `ToolResult`，命中表路径在 `artifacts`；解析函数返回含位置、链方向和分数的命中字典。
- **示例代码：**

```python
from dnakit.bioinformatics.alignment import build_profile, search_profile, parse_nhmmer_table

model = build_profile("aligned.sto", "new_profile", alphabet="dna")
result = search_profile(model.artifacts["profile.hmm"], "targets.fa", "new_hmm_search")
for hit in parse_nhmmer_table(result.artifacts["hits.tbl"]):
    print(hit["target"], hit["start"], hit["end"], hit["strand"])
```

- **示例结果：**

命中为字典，含 `target`、`query`、`start/end`、`strand`、`evalue`、`score`；坐标为 0-based 半开区间。外部程序示例未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#alignment)。

<!-- BEGIN GENERATED API: blast -->

<span id="blast"></span>

## 9) 使用远程 BLAST 检索相似蛋白质序列 {#dnakit.databases.blast_protein_search}

- **作用：** 使用远程 BLAST 检索相似蛋白质序列。
- **处理方法：** 通过 `NCBIBlast` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.blast_protein_search(*, sequence, database=None, expect=None, hitlist_size=None, verbose=True)`
- **输入：** 必须提供 `sequence`。
- **输出：** 返回字典，包含 `query_id`, `query_length`, `database`, `algorithm` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.blast_protein_search(
    sequence='MKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL',
    database='swissprot',
    hitlist_size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "query_id": "<string>",
  "query_length": "<integer>",
  "database": "<string>",
  "algorithm": "<string>",
  "alignments": [
    {
      "hit_id": "<string>",
      "hit_def": "<string>",
      "hit_length": "<integer>",
      "hsps": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/09_search.md#dnakit.databases.blast_protein_search-parameters)

<!-- END GENERATED API: blast -->
