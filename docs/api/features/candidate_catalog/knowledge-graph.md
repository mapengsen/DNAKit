# 知识图谱与数据溯源

本类收录 **18** 项有文档的功能。其中 17 项处理以普通三元组/边列表表示的图（纯 Python，不使用图数据库或网络服务）：这 17 项全部在本地计算出结果，没有一项返回 `not_implemented`。剩下的 `pseudonymize_genomic_metadata_fields` 与图无关：它对元数据映射中的敏感字段做假名化（基于键名和模式，并非去标识化），见其自己的小节。图输入可以是映射、三元组列表、带分隔符的文本或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径，因此所有函数既能在 Python 中调用，也能通过 JSON/MCP 调用。请注意 `apply_access_policy` 的默认行为：没有密级标签的项会被移除，除非传入 `unclassified='allow'`。下面的示例中，较长的数组和文本用 `...` 缩短。

## 1) add_graph_relationship

- **作用：** 向知识图谱添加一条或多条关系（边），返回新图；输入不会被修改。

- **处理方法：** 加载图（输入为 `None` 时从空图开始），添加给定的边及其属性，除非 `create_missing_nodes` 为假，否则自动创建缺失节点。边由 (source, relation, target) 标识；重复边按 `on_duplicate` 处理。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.add_graph_relationship(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）；或 `None`。`options`：`subject`、`predicate`、`object` 和可选的 `properties`（单条边），或 `relationships`（三元组/边映射列表）；`create_missing_nodes`（`True`；为 `False` 时未知节点会报错）、`on_duplicate`（默认 `ignore`，另有 `error`、`update`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`added`（新增的边）、`already_present`（已存在的边）和 `graph`（`directed`、`nodes`、`edges`）。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import add_graph_relationship

graph = [["BRCA1", "interacts_with", "TP53"]]
result = add_graph_relationship(graph, options={"subject": "TP53", "predicate": "regulates",
                                                "object": "CDKN1A", "properties": {"source": "demo"}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "add_graph_relationship",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 37,
    "sequence_bases": 19,
    "lengths": [19],
    "sha256": "85b97b4918132e5d079845aba613f064b46e4ce8b6976c331fd4db8d48a2ac83"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 3,
    "n_edges": 2,
    "relation_counts": {"interacts_with": 1, "regulates": 1},
    "added": [
      {"source": "TP53", "relation": "regulates", "target": "CDKN1A"}
    ],
    "already_present": [],
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "BRCA1"},
        {"id": "TP53"},
        {"id": "CDKN1A"}
      ],
      "edges": [
        {"source": "BRCA1", "relation": "interacts_with", "target": "TP53"},
        {"source": "TP53", "relation": "regulates", "target": "CDKN1A"}
      ]
    }
  }
}
```

## 2) aggregate

- **作用：** 把知识图谱汇总为统计量：关系计数、度统计、密度、连通分量、节点类型和度最高的枢纽节点。

- **处理方法：** 度 = 图中边集合里的入边 + 出边（自环计两次）；密度 = 边数 / 可能的边数（不含自环）；连通分量为弱连通分量；`top_hubs` 为度最高的 `top_k` 个节点；节点类型按各节点的 `type_property` 计数（缺失记为 `unknown`）。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.aggregate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`top_k`（10）、`type_property`（`type`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`degree`（`total`、`in`、`out`，各含最小/最大/均值/中位数）、`density`、`n_isolated_nodes`、`n_self_loops`、`n_weakly_connected_components`、`largest_component_size`、`node_type_counts`、`top_hubs` 和 `degree_note`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import aggregate

triples = [["BRCA1", "associated_with", "breast_cancer"],
           ["TP53", "associated_with", "breast_cancer"],
           ["BRCA1", "interacts_with", "TP53"],
           ["TP53", "regulates", "CDKN1A"]]
result = aggregate(triples, options={"top_k": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "aggregate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 163,
    "sequence_bases": 93,
    "lengths": [30, 28, 19, 16],
    "sha256": "489d99ad8c84844d7ffbbbfc85328cced425e1674a460d69d9f49979ac6043b0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 4,
    "n_edges": 4,
    "relation_counts": {"associated_with": 2, "interacts_with": 1, "regulates": 1},
    "degree": {
      "total": {"min": 1, "max": 3, "mean": 2.0, "median": 2.0},
      "in": {"min": 0, "max": 2, "mean": 1.0, "median": 1.0},
      "out": {"min": 0, "max": 2, "mean": 1.0, "median": 1.0}
    },
    "density": 0.3333333333333333,
    "n_isolated_nodes": 0,
    "n_self_loops": 0,
    "n_weakly_connected_components": 1,
    "largest_component_size": 4,
    "node_type_counts": {"unknown": 4},
    "top_hubs": [
      {"id": "TP53", "degree": 3, "in_degree": 1, "out_degree": 2},
      {"id": "BRCA1", "degree": 2, "in_degree": 0, "out_degree": 2}
    ],
    "degree_note": "degree = in + out edges of the graph's edge set; a self-loop counts twice"
  }
}
```

## 3) apply_access_policy

- **作用：** 对知识图谱应用按密级的访问策略：移除密级高于许可级别的节点和边，并脱敏属性。

- **处理方法：** 每个节点和边在属性 `classification_property` 中带一个密级标签，`levels` 按由低到高排列这些标签。密级不高于 `clearance` 的项被保留；端点节点已被移除的边也被移除，随后从保留项中删去 `redact_properties` 列出的属性。**需要注意的默认值：**`clearance` 默认为最低级别（`public`），`unclassified` 默认为 `deny`，即没有密级标签的节点或边一律被移除。因此，没有任何密级标签的图在默认选项下会返回空图。要保留未分级的项，请传 `unclassified='allow'`。不在 `levels` 中的标签会报 `UNKNOWN_CLASSIFICATION`。它只是过滤数据，不是认证系统。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.apply_access_policy(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`clearance`（默认最低级别）、`levels`（由低到高，默认 `['public', 'internal', 'restricted', 'secret']`）、`classification_property`（`classification`，从节点和边的属性读取）、`unclassified`（默认 `deny`：没有密级的项被移除；或 `allow`）、`deny_nodes`（无论如何都移除的 id）、`redact_properties`（从保留项中删去的属性名）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含过滤后图的 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`clearance`、`levels`、`n_nodes_removed`、`removed_nodes`（至多 1000 个）、`n_edges_removed`、说明未分级项被移除还是保留的 `policy` 句子，以及 `graph`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import apply_access_policy

graph = {
    "nodes": [
        {"id": "gene1", "classification": "public"},
        {"id": "patient1", "classification": "restricted", "name": "A. Person"},
        {"id": "sample1", "classification": "internal", "note": "freezer 3"},
        {"id": "gene2"},                                  # no classification
    ],
    "edges": [
        {"source": "gene1", "relation": "expressed_in", "target": "sample1", "classification": "internal"},
        {"source": "patient1", "relation": "donated", "target": "sample1", "classification": "restricted"},
        {"source": "gene2", "relation": "expressed_in", "target": "sample1", "classification": "public"},
    ],
}
# default unclassified="deny": the unclassified node gene2 (and its edge) is removed;
# unclassified="allow" would keep it
result = apply_access_policy(graph, options={"clearance": "internal", "redact_properties": ["note"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "apply_access_policy",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 525,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "bddd261b368cafdb02d71b3e81fb54de312f101516629a16bca4a4689700b50c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 2,
    "n_edges": 1,
    "relation_counts": {"expressed_in": 1},
    "clearance": "internal",
    "levels": ["public", "internal", "restricted", "secret"],
    "n_nodes_removed": 2,
    "removed_nodes": ["patient1", "gene2"],
    "n_edges_removed": 2,
    "policy": "keep items with classification <= internal; unclassified items are removed",
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "gene1", "classification": "public"},
        {"id": "sample1", "classification": "internal"}
      ],
      "edges": [
        {"source": "gene1", "relation": "expressed_in", "target": "sample1", "classification": "internal"}
      ]
    }
  }
}
```

## 4) build_data_lineage

- **作用：** 根据溯源记录构建数据血缘图（输出 derived_from 输入），并报告根节点、叶节点和拓扑顺序。

- **处理方法：** 每条记录添加从输出指向各输入的边，关系为 `derived_from`，步骤名存在边上。根节点是没有输入的节点，叶节点是没有其他节点由它派生的节点；拓扑顺序把输入排在输出之前。存在环时会报错。给出 `trace` 时，返回该节点的全部上游祖先及其距离（广度优先）。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.build_data_lineage(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：记录列表 `{'output': id（或 'target' / 'dataset' / 'id'）, 'inputs': [ids]（或 'sources' / 'derived_from'）, 'step': 名称（或 'tool' / 'process'）}`。`options`：`trace`（节点 id）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`roots`、`leaves`、`topological_order`、`edge_meaning`、`graph`，给出 `trace` 时另有 `upstream_of`（`node`、含 `id` 与 `distance` 的 `ancestors`）。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import build_data_lineage

records = [
    {"output": "aligned.bam", "inputs": ["reads.fastq", "ref.fa"], "step": "bwa"},
    {"output": "variants.vcf", "inputs": ["aligned.bam"], "step": "caller"},
]
result = build_data_lineage(records, options={"trace": "variants.vcf"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "build_data_lineage",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 152,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "59c5f344703985c993aec3978b80fbd06ebe1baf7a90c3085ec01f97dda6336e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 4,
    "n_edges": 3,
    "relation_counts": {"derived_from": 3},
    "roots": ["reads.fastq", "ref.fa"],
    "leaves": ["variants.vcf"],
    "topological_order": ["reads.fastq", "ref.fa", "aligned.bam", "variants.vcf"],
    "edge_meaning": "source was derived_from target (output -> input)",
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "aligned.bam"},
        {"id": "reads.fastq"},
        {"id": "ref.fa"},
        {"id": "variants.vcf"}
      ],
      "edges": [
        {"source": "aligned.bam", "relation": "derived_from", "target": "reads.fastq", "step": "bwa"},
        {"source": "aligned.bam", "relation": "derived_from", "target": "ref.fa", "step": "bwa"},
        {"source": "variants.vcf", "relation": "derived_from", "target": "aligned.bam", "step": "caller"}
      ]
    },
    "upstream_of": {
      "node": "variants.vcf",
      "ancestors": [
        {"id": "aligned.bam", "distance": 1},
        {"id": "reads.fastq", "distance": 2},
        {"id": "ref.fa", "distance": 2}
      ]
    }
  }
}
```

## 5) build_knowledge_graph

- **作用：** 由三元组或边列表构建知识图谱（节点唯一，重复边合并）。

- **处理方法：** 仅在边中出现的节点会被创建；重复的 (source, relation, target) 边被合并，重复次数存入边属性 `count`。除非 `drop_self_loops` 为真，否则保留自环。无向图（`directed=False`）中，边的两端不区分顺序。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.build_knowledge_graph(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`duplicate_nodes`、`duplicate_edges`、`self_loops`、`dropped_self_loops` 和 `graph`（`directed`、`nodes`、`edges`）。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import build_knowledge_graph

triples = [["A", "binds", "B"], ["A", "binds", "B"], ["B", "binds", "B"], ["C", "binds", "A"]]
result = build_knowledge_graph(triples)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "build_knowledge_graph",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 84,
    "sequence_bases": 28,
    "lengths": [7, 7, 7, 7],
    "sha256": "c25bb514c2c7897ece0f8c5a3d04e6e45b3f6c3365135ef6c14d7ce908fc7730"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 3,
    "n_edges": 3,
    "relation_counts": {"binds": 3},
    "duplicate_nodes": 0,
    "duplicate_edges": 1,
    "self_loops": 1,
    "dropped_self_loops": 0,
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "A"},
        {"id": "B"},
        {"id": "C"}
      ],
      "edges": [
        {"source": "A", "relation": "binds", "target": "B", "count": 2},
        {"source": "B", "relation": "binds", "target": "B"},
        {"source": "C", "relation": "binds", "target": "A"}
      ]
    }
  }
}
```

## 6) compare_database_versions

- **作用：** 比较知识图谱的两个版本：新增、删除和变更的节点与边。

- **处理方法：** 节点以 id 标识，边以 (source, relation, target) 标识；属性不同的节点或边视为“变更”。同时给出两个基于内容的版本号（见 `version`），二者相等时 `identical` 为真。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.compare_database_versions(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'old': 图, 'new': 图}`、由图映射或文件组成的二元组 `[old, new]`，或一个图加 `options['other']` 作为新版本。图的形式同 `build_knowledge_graph`。`options`：`max_items`（1000，限制每个列表的长度）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `identical`、`old_version_id`、`new_version_id`、`counts`（节点/边的新增、删除、变更数）、`added_nodes`、`removed_nodes`、`changed_nodes`（`id`、`old`、`new` 属性）、`added_edges`、`removed_edges`、`changed_edges` 和 `truncated`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import compare_database_versions

old = {"nodes": [{"id": "A", "score": 1}, "B"], "edges": [["A", "binds", "B"]]}
new = {"nodes": [{"id": "A", "score": 2}, "B", "C"], "edges": [["A", "binds", "B"], ["B", "binds", "C"]]}
result = compare_database_versions({"old": old, "new": new})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "compare_database_versions",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 190,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7043195422a8b0e174b381b6c9cd45951da6b947604adc976bf3340e3b3ace40"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "identical": false,
    "old_version_id": "sha256:adc9d9cb04bcc4a5d1360dfe54789cad53df0e170b695aa43207f0213fc85465",
    "new_version_id": "sha256:a2968fe2e09bbc3df6b021824c3660ab4deac76c0899257c4956f3b33c41c06b",
    "counts": {
      "nodes_added": 1,
      "nodes_removed": 0,
      "nodes_changed": 1,
      "edges_added": 1,
      "edges_removed": 0,
      "edges_changed": 0
    },
    "added_nodes": ["C"],
    "removed_nodes": [],
    "changed_nodes": [
      {
        "id": "A",
        "old": {"score": 1},
        "new": {"score": 2}
      }
    ],
    "added_edges": [
      {"source": "B", "relation": "binds", "target": "C"}
    ],
    "removed_edges": [],
    "changed_edges": [],
    "truncated": false
  }
}
```

## 7) create

- **作用：** 创建新的知识图谱，可选地预置节点和边。

- **处理方法：** 从空图（输入为 `None`）或 `build_knowledge_graph` 接受的任意图输入开始，再加入 `nodes` 与 `edges` 选项，并把 `name`、`description` 存为图的元数据。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.create(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`None`（空图），或 `build_knowledge_graph` 接受的任意图输入。`options`：`name`、`description`、`directed`（`True`）、`nodes`（id 列表或 `{'id': ..., 属性}`）、`edges`（三元组或边映射列表）、`default_relation`。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`duplicate_nodes`、`duplicate_edges`、`self_loops`、`dropped_self_loops` 和 `graph`（给出时含 `name` 与 `description`）。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import create

result = create(None, options={"name": "demo", "nodes": ["A", {"id": "B", "type": "gene"}],
                               "edges": [["A", "binds", "B"]]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "create",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 2,
    "n_edges": 1,
    "relation_counts": {"binds": 1},
    "duplicate_nodes": 0,
    "duplicate_edges": 0,
    "self_loops": 0,
    "dropped_self_loops": 0,
    "graph": {
      "name": "demo",
      "directed": true,
      "nodes": [
        {"id": "A"},
        {"id": "B", "type": "gene"}
      ],
      "edges": [
        {"source": "A", "relation": "binds", "target": "B"}
      ]
    }
  }
}
```

## 8) pseudonymize_genomic_metadata_fields

- **作用：** 对元数据映射（或映射列表）做字段假名化：在任意嵌套深度替换、删除或用带密钥的哈希处理识别性键的值，并替换文本中的电子邮件地址、电话号码和类似社保号的模式。这是基于键名和模式的假名化，并不保证完成去标识化。它取代了早先的一个函数——那个函数只用无盐 SHA-256 处理六个顶层键，名称承诺的比实际做到的多，因此已被移除。

- **处理方法：** 在 `src/dnakit/_metadata_pseudonymize.py` 中只实现一次，并与 `dnakit.workflow_security.pseudonymize_metadata_fields` 共用。输入不会被修改，对深拷贝做转换。当一个键的名称转为小写并去掉字母和数字以外的所有字符后（因此 `Patient_ID`、`patient id`、`PATIENTID` 视为相同）出现在默认列表中时，该键被视为敏感键：`name`、`first/last/full/given/family/middle name`、`surname`、`patient`、`patient id`、`patient name`、`mrn`、`medical record number`、`ssn`、`social security number`、`email`、`email address`、`phone`、`phone number`、`telephone`、`mobile`、`fax`、`address`、`home/street address`、`street`、`city`、`zip`、`zip code`、`postal code`、`postcode`、`dob`、`date of birth`、`birth date`、`birthday`、`national id`、`passport`、`passport number`、`insurance`、`insurance id/number`、`account number`、`ip address`、`device id`、`photo`、`photograph`、`biometric(s)`、`license number`、`driver license`。匹配是对规范化名称的精确匹配（`sample_name` 不会命中）；`options['keys']` 增加名称，`options['keep']` 去掉名称（keep 优先）。敏感键的值整体被替换，值为映射或列表时也一样。`mode`：`redact`（默认）写入 `"[REDACTED]"`，不需要任何密钥；`drop` 删除该键；`hash` 写入以 `options['salt']`（必填，至少 8 个字符的字符串，否则抛出 `ConfigurationError`）为密钥、对 UTF-8 值（字符串按原样，其他值用紧凑的、键排序的 JSON）计算的 HMAC-SHA256 的前 `digest_chars`（默认 16，范围 8 到 64）个十六进制字符。同一盐值下哈希是确定的，所以相同的值仍可关联，知道盐值的人可以猜出低熵的值；请保密盐值。启用 `scan_text`（默认 `True`）时，任何位置的字符串值（包括保留键下的值和列表中的值）都会被扫描电子邮件地址、类似社保号的数字（`123-45-6789`）、国际电话号码（`+44 20 7946 0958`）和北美格式的 10 位号码；命中部分变为 `[REDACTED:email]`、`[REDACTED:ssn]`、`[REDACTED:phone]`。电话模式也会命中单纯的 10 位编号，因此可能过度脱敏。它不做的事：不对日期做平移或泛化，不泛化超过 89 岁的年龄和小地理单位，不识别自由文本中的姓名和其他准标识符，不改写本身含有标识符的字典键，也不对照 HIPAA 或 GDPR 的标识符清单检查。共享前请人工复核输出。datetime 会转为 ISO 字符串；其他非 JSON 类型（集合、bytes、对象）、超过 64 层的嵌套或循环引用会抛出 `InputFormatError`。

- **API：** `dnakit.knowledge_graph.pseudonymize_genomic_metadata_fields(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：映射、映射列表，或其中之一的 JSON 文本；其他输入抛出 `InputFormatError`（错误码 `INVALID_METADATA_INPUT`、`INVALID_METADATA_JSON`、`UNSUPPORTED_METADATA_VALUE`、`METADATA_TOO_DEEP`，均带 `hint`）。`options`：`mode`（`"redact"` | `"drop"` | `"hash"`，默认 `"redact"`）、`salt`（至少 8 个字符的字符串，`hash` 必填）、`digest_chars`（8 到 64，默认 16，仅 `hash`）、`keys`（额外的敏感键名，字符串或列表）、`keep`（永不视为敏感的键名）、`scan_text`（默认 `True`）、`backend`（仅作标签）。无效的选项抛出 `ConfigurationError`（`PSEUDONYMIZE_SALT_REQUIRED`、`INVALID_PSEUDONYMIZE_MODE`、`INVALID_DIGEST_CHARS`、`INVALID_PSEUDONYMIZE_KEYS`、`INVALID_PSEUDONYMIZE_OPTIONS`）。

- **输出：** 标准封装：`status`（`"ok"`）、`module`（`"knowledge_graph"`）、`operation`（`"pseudonymize_genomic_metadata_fields"`）、`provenance`（`implementation` 为 `"native"`，另有 `backend`）、`input_summary`、`warnings` 和 `result`。`input_summary` 特意只含 `record_count` 和 `n_keys`（访问过的映射键个数）：不含长度，也不含未受保护输入的 SHA-256。`result` 含 `metadata`（转换后的映射或列表；JSON 文本输入以解析后的结构返回）、`mode`、`n_redacted`（= `n_key_matches` + `n_text_matches`）、`n_key_matches`、`n_text_matches`、`text_matches_by_pattern`、`redacted_paths`（类 JSON 路径的键路径，如 `$.visits[0].dob`，绝不含值；最多 10000 条，`redacted_paths_truncated` 表示是否还有更多）和 `warnings`。警告总是声明：该方法基于键名和模式，并不保证完成去标识化，不对日期、89 岁以上年龄或小地理单位做平移或泛化，自由文本仍可能识别出个人，共享前必须复核输出；同一份列表同时出现在 `result.warnings` 和封装的 `warnings` 中。输出的任何位置都不会出现盐值或原始值。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import pseudonymize_genomic_metadata_fields

metadata = {
    "name": "Ada Lovelace",
    "patient_id": "P-001",
    "Email": "ada@example.org",
    "age": 36,
    "visit": {
        "name": "Ada Lovelace",
        "note": "call +1 (415) 555-0100 or ada@example.org, SSN 123-45-6789",
    },
}

# default: values of sensitive keys at any depth become "[REDACTED]"; patterns inside text are replaced
result = pseudonymize_genomic_metadata_fields(metadata)
print(json.dumps(result, indent=2))

# keyed hash: deterministic tokens that depend on a secret salt (required, >= 8 characters)
hashed = pseudonymize_genomic_metadata_fields(metadata, options={"mode": "hash", "salt": "example-salt-123", "digest_chars": 12})
print(json.dumps(hashed["result"]["metadata"], indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "pseudonymize_genomic_metadata_fields",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "n_keys": 7
  },
  "warnings": [
    "Pseudonymisation is key-name and pattern based only. It is NOT a guarantee of de-identification and does not by itself satisfy HIPAA or GDPR requirements.",
    "Dates are not shifted or generalised; ages over 89 and small geographic units are not generalised.",
    "Free text may still identify a person: only e-mail, phone and SSN-like patterns are scanned, names and other identifiers inside text are not detected; dictionary keys that themselves contain identifiers are not rewritten.",
    "Review the output manually before sharing it."
  ],
  "status": "ok",
  "result": {
    "metadata": {
      "name": "[REDACTED]",
      "patient_id": "[REDACTED]",
      "Email": "[REDACTED]",
      "age": 36,
      "visit": {
        "name": "[REDACTED]",
        "note": "call [REDACTED:phone] or [REDACTED:email], SSN [REDACTED:ssn]"
      }
    },
    "mode": "redact",
    "n_redacted": 7,
    "n_key_matches": 4,
    "n_text_matches": 3,
    "text_matches_by_pattern": {
      "email": 1,
      "phone": 1,
      "ssn": 1
    },
    "redacted_paths": [
      "$.name",
      "$.patient_id",
      "$.Email",
      "$.visit.name",
      "$.visit.note"
    ],
    "redacted_paths_truncated": false,
    "warnings": [
      "Pseudonymisation is key-name and pattern based only. It is NOT a guarantee of de-identification and does not by itself satisfy HIPAA or GDPR requirements.",
      "Dates are not shifted or generalised; ages over 89 and small geographic units are not generalised.",
      "Free text may still identify a person: only e-mail, phone and SSN-like patterns are scanned, names and other identifiers inside text are not detected; dictionary keys that themselves contain identifiers are not rewritten.",
      "Review the output manually before sharing it."
    ]
  }
}
{
  "name": "69272bf887f1",
  "patient_id": "61df4f560250",
  "Email": "0e934605ef5f",
  "age": 36,
  "visit": {
    "name": "69272bf887f1",
    "note": "call [REDACTED:phone] or [REDACTED:email], SSN [REDACTED:ssn]"
  }
}
```

## 9) export

- **作用：** 把知识图谱序列化为 JSON、JSONL、TSV（仅边三元组）或 GraphML 文本，可选地写入文件。

- **处理方法：** GraphML 把属性值存为字符串；TSV 只保存边三元组，因此节点和边的属性以及孤立节点都会丢失（`result` 中有 `note` 提示）。文件以原子方式写入，除非 `overwrite` 为真，否则不会覆盖已有文件。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`format`（`json`、`jsonl`、`tsv`、`graphml`；默认由 `output` 后缀决定，否则 `json`）、`output`（路径；不给则文本以 `content` 返回）、`overwrite`（`False`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`format`、`bytes`、`content`（未给出 `output` 时）、`output_path`，TSV 另有 `note`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import export

result = export([["A", "binds", "B"], ["B", "regulates", "C"]], options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "export",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 46,
    "sequence_bases": 18,
    "lengths": [7, 11],
    "sha256": "212a4cf1305741b038617d44f25ccc4ae03c9d9744a3126339d18c125bd523e1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 3,
    "n_edges": 2,
    "relation_counts": {"binds": 1, "regulates": 1},
    "format": "tsv",
    "bytes": 47,
    "sha256": "933baf6d7d3b0b7f31989f777d88bf0710f0dd6cd1f33d33bd66e17622d95b86",
    "content": "source\trelation\ttarget\nA\tbinds\tB\nB\tregulates\tC\n",
    "output_path": null,
    "note": "TSV holds edge triples only: node/edge properties are not written."
  }
}
```

## 10) link

- **作用：** 关联实体：在明确给出的节点对之间建边，或在某个节点属性取值相同的所有节点之间建边。

- **处理方法：** 给出 `pairs` 时，对每个 [a, b] 添加关系为 `relation` 的边（缺失节点会被创建）。给出 `by` 时，该属性取值相等（忽略大小写）的节点两两相连，每组最多 `max_group_size` 个节点。已存在的边不会重复添加。`pairs` 与 `by` 必须恰好给出一个。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.link(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`relation`（`linked_to`）；`pairs`（`[a, b]` 列表）与 `by`（节点属性名）二选一；`max_group_size`（100）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`n_created`、`created`（新增的边）和 `graph`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import link

graph = {"nodes": [{"id": "g1", "pathway": "p53"}, {"id": "g2", "pathway": "P53"}, {"id": "g3", "pathway": "wnt"}],
         "edges": []}
result = link(graph, options={"by": "pathway", "relation": "same_pathway"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "link",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 120,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "2e0b94eb2020347dda8671110b39634d9a155408699c6d912989561bc9b331b3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 3,
    "n_edges": 1,
    "relation_counts": {"same_pathway": 1},
    "n_created": 1,
    "created": [
      {"source": "g1", "relation": "same_pathway", "target": "g2"}
    ],
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "g1", "pathway": "p53"},
        {"id": "g2", "pathway": "P53"},
        {"id": "g3", "pathway": "wnt"}
      ],
      "edges": [
        {"source": "g1", "relation": "same_pathway", "target": "g2"}
      ]
    }
  }
}
```

## 11) map

- **作用：** 通过查找表映射（重命名）节点标识符和关系名，并合并映射到同一标识符的节点。

- **处理方法：** 节点 id 通过 `mapping` 替换，关系名通过 `relation_mapping` 替换。映射到同一新 id 的节点被合并（合并后的节点保留第一个节点的属性），变得相同的边也被合并并计入 `count`。`mapping` 中没有的节点按 `unmapped` 保留、丢弃或报错。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.map(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`mapping`（旧 id -> 新 id 的字典、行 `{'from', 'to'}` / `[from, to]`，或两列 TSV/CSV 路径）、`relation_mapping`（形式相同，用于关系名）、`unmapped`（默认 `keep`，另有 `drop`、`error`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`n_ids_renamed`、`n_nodes_merged`、`n_nodes_dropped` 和 `graph`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import map as map_graph

graph = [["ENSG1", "interacts_with", "ENSG2"], ["ENSG3", "interacts_with", "ENSG2"]]
result = map_graph(graph, options={"mapping": {"ENSG1": "BRCA1", "ENSG2": "TP53", "ENSG3": "BRCA1"}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "map",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 76,
    "sequence_bases": 42,
    "lengths": [21, 21],
    "sha256": "d2b0ea1fb95e0ba2ea8201455878623be0a40403eeb222c33c7733cd99c36096"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 2,
    "n_edges": 1,
    "relation_counts": {"interacts_with": 1},
    "n_ids_renamed": 3,
    "n_nodes_merged": 1,
    "n_nodes_dropped": 0,
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "BRCA1"},
        {"id": "TP53"}
      ],
      "edges": [
        {"source": "BRCA1", "relation": "interacts_with", "target": "TP53", "count": 2}
      ]
    }
  }
}
```

## 12) query_graph_neighbors

- **作用：** 用广度优先搜索查询节点的邻域，可选地给出到目标节点的最短路径和诱导子图。

- **处理方法：** 从 `node` 出发做深度至多为 `depth` 的广度优先搜索，按 `direction` 沿边前进（无向图忽略），若给出 `relations` 则只走这些关系；每个邻居给出距离、一条最短路径和到达它所经的边。给出 `target` 时，返回（同样方向和关系规则下）边数最少的最短路径，长度至多 `max_path_length`。`node` 未知时会报错。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.query_graph_neighbors(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`node`（必需）、`depth`（1，最大 20）、`direction`（默认 `both`，另有 `out`、`in`）、`relations`（仅限这些关系名）、`limit`（1000）、`target`、`max_path_length`（50）、`include_subgraph`（`False`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `node`、`depth`、`direction`、`n_neighbors`、`neighbors`（`id`、`distance`、`path`、`via_edge`）、`truncated`；给出 `target` 时另有 `shortest_path`（`found`、`length`、`nodes`、`edges`），给出 `include_subgraph` 时另有诱导子图 `subgraph`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import query_graph_neighbors

triples = [["BRCA1", "associated_with", "breast_cancer"],
           ["TP53", "associated_with", "breast_cancer"],
           ["BRCA1", "interacts_with", "TP53"],
           ["TP53", "regulates", "CDKN1A"]]
result = query_graph_neighbors(triples, options={"node": "BRCA1", "depth": 2,
                                                  "direction": "out", "target": "CDKN1A"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "query_graph_neighbors",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 163,
    "sequence_bases": 93,
    "lengths": [30, 28, 19, 16],
    "sha256": "489d99ad8c84844d7ffbbbfc85328cced425e1674a460d69d9f49979ac6043b0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "node": "BRCA1",
    "depth": 2,
    "direction": "out",
    "n_neighbors": 3,
    "neighbors": [
      {
        "id": "breast_cancer",
        "distance": 1,
        "path": ["BRCA1", "breast_cancer"],
        "via_edge": {"source": "BRCA1", "relation": "associated_with", "target": "breast_cancer"}
      },
      {
        "id": "TP53",
        "distance": 1,
        "path": ["BRCA1", "TP53"],
        "via_edge": {"source": "BRCA1", "relation": "interacts_with", "target": "TP53"}
      },
      {
        "id": "CDKN1A",
        "distance": 2,
        "path": ["BRCA1", "TP53", "CDKN1A"],
        "via_edge": {"source": "TP53", "relation": "regulates", "target": "CDKN1A"}
      }
    ],
    "truncated": false,
    "shortest_path": {
      "found": true,
      "length": 2,
      "nodes": ["BRCA1", "TP53", "CDKN1A"],
      "edges": [
        {"source": "BRCA1", "relation": "interacts_with", "target": "TP53"},
        {"source": "TP53", "relation": "regulates", "target": "CDKN1A"}
      ]
    }
  }
}
```

## 13) read

- **作用：** 从文件或文本（JSON、JSONL、TSV/CSV 三元组、GraphML）读取知识图谱，并连同摘要一起返回。

- **处理方法：** 按图输入的各种形式解析，并像 `build_knowledge_graph` 一样构建图。含 DOCTYPE/ENTITY 声明的 GraphML 会被拒绝，大于 100 MB 的文件也会被拒绝。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`duplicate_nodes`、`duplicate_edges`、`self_loops`、`dropped_self_loops` 和 `graph`。

- **示例代码：**

```python
import json
from pathlib import Path
from dnakit.knowledge_graph import read

Path("edges.tsv").write_text("BRCA1\tinteracts_with\tTP53\nTP53\tregulates\tCDKN1A\n")
result = read("edges.tsv")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "read",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 9,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "91d122db62edbb273e5b0672dcbc62a8d0d0f4470a6ea420f37fc8028783bfbb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 3,
    "n_edges": 2,
    "relation_counts": {"interacts_with": 1, "regulates": 1},
    "duplicate_nodes": 0,
    "duplicate_edges": 0,
    "self_loops": 0,
    "dropped_self_loops": 0,
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "BRCA1"},
        {"id": "TP53"},
        {"id": "CDKN1A"}
      ],
      "edges": [
        {"source": "BRCA1", "relation": "interacts_with", "target": "TP53"},
        {"source": "TP53", "relation": "regulates", "target": "CDKN1A"}
      ]
    }
  }
}
```

## 14) search

- **作用：** 按节点 id 及字符串/数值属性中的文本搜索知识图谱的节点。

- **处理方法：** 把查询词与节点 id 及每个字符串/数值属性（若给出 `fields` 则仅限这些字段）按 `match`（子串、完全匹配或前缀；除非 `case_sensitive`，否则不区分大小写）比较。可把节点限制为某个 `node_type`。匹配结果保持图中的节点顺序，并给出命中的字段。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.search(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`query`（必需）、`match`（默认 `contains`，另有 `exact`、`prefix`）、`case_sensitive`（`False`）、`fields`（要搜索的属性名，包含 `id`）、`node_type` 与 `type_property`（`type`）、`limit`（100）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `query`、`match`、`n_matches`、`matches`（节点 `id`、`matched_fields` 与节点属性）和 `truncated`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import search

graph = {"nodes": [{"id": "BRCA1", "type": "gene", "description": "DNA repair"},
                   {"id": "TP53", "type": "gene", "description": "tumour suppressor"},
                   {"id": "repair_pathway", "type": "pathway"}], "edges": []}
result = search(graph, options={"query": "repair"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "search",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 199,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "513a79306bb4d244ac0b7536b5cefdf7948dbef629b50172481e4083083cb3f4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "query": "repair",
    "match": "contains",
    "n_matches": 2,
    "matches": [
      {"id": "BRCA1", "matched_fields": ["description"], "type": "gene", "description": "DNA repair"},
      {"id": "repair_pathway", "matched_fields": ["id"], "type": "pathway"}
    ],
    "truncated": false
  }
}
```

## 15) snapshot_database_version

- **作用：** 为知识图谱创建快照（规范化内容加版本号），可选地保存为 JSON 文件。

- **处理方法：** 快照由规范化（排序后）的图及其基于内容的 `version_id`（同 `version`）、可选的 `label` 以及（仅在给出时才有的）`timestamp` 组成，因此同一个图的多次快照逐字节相同。文件以原子方式写入，除非 `overwrite` 为真，否则不会覆盖已有文件。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.snapshot_database_version(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`label`、`timestamp`（仅在给出时保存）、`output`（路径；不给则返回快照）、`overwrite`（`False`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `version_id`、`directed`、`n_nodes`、`n_edges`、`relation_counts`、`snapshot`（`format`、`format_version`、`version_id`、`label`、`graph`）、`bytes`、`content`（JSON 文本，未给出 `output` 时）和 `output_path`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import snapshot_database_version

result = snapshot_database_version([["A", "binds", "B"]], options={"label": "v1"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "snapshot_database_version",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 21,
    "sequence_bases": 7,
    "lengths": [7],
    "sha256": "f7d05ee5ab5f2f13c73b45300a8c5cf4a16555e0732720d87be5e845280e9713"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "version_id": "sha256:14f2eff389d5a427f97439a4bd834d26d9351a6185dea4968884bd17f501baf3",
    "directed": true,
    "n_nodes": 2,
    "n_edges": 1,
    "relation_counts": {"binds": 1},
    "snapshot": {
      "format": "dnakit.knowledge_graph.snapshot",
      "format_version": 1,
      "version_id": "sha256:14f2eff389d5a427f97439a4bd834d26d9351a6185dea4968884bd17f501baf3",
      "label": "v1",
      "graph": {
        "directed": true,
        "nodes": [
          {"id": "A"},
          {"id": "B"}
        ],
        "edges": [
          {"source": "A", "relation": "binds", "target": "B"}
        ]
      }
    },
    "bytes": 423,
    "sha256": "e646a6ac4e504dbb35cf96c57df68e48f91687b4901645125248f1c6148fe095",
    "content": "{\n  \"format\": \"dnakit.knowledge_graph.snapshot\",\n  \"format_version\": 1,\n  \"version_id\": \"sha256:14f2eff389d5a427f97439a4bd834d26d9351a6185dea4968884bd17f501baf3...",
    "output_path": null
  }
}
```

## 16) update

- **作用：** 对知识图谱应用一组变更：删除边或节点（连同其边），并设置节点或边的属性。

- **处理方法：** 返回应用了所请求变更的新图；输入不会被修改。删除节点会同时删除它的所有边。除非 `ignore_missing` 为真，未知 id 会报错。至少需要给出一个变更选项。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.update(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`（至少一个）：`remove_edges`（`[source, relation, target]` 列表）、`remove_nodes`（id 列表）、`set_node_properties`（`{id: {属性: 值}}`）、`set_edge_properties`（`{source, relation, target, properties}` 列表）、`ignore_missing`（`False`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `directed`、`n_nodes`、`n_edges` 与 `relation_counts`、`edges_removed`、`nodes_removed`、`node_properties_set`、`edge_properties_set` 和 `graph`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import update

graph = [["A", "binds", "B"], ["B", "binds", "C"]]
result = update(graph, options={"remove_nodes": ["C"], "set_node_properties": {"A": {"type": "gene"}}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "update",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 42,
    "sequence_bases": 14,
    "lengths": [7, 7],
    "sha256": "80f96dc5cc4579423529eb7d64469b6493b6959097d5721000531a1bc7f118c0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "directed": true,
    "n_nodes": 2,
    "n_edges": 1,
    "relation_counts": {"binds": 1},
    "edges_removed": 1,
    "nodes_removed": 1,
    "node_properties_set": 1,
    "edge_properties_set": 0,
    "graph": {
      "directed": true,
      "nodes": [
        {"id": "A", "type": "gene"},
        {"id": "B"}
      ],
      "edges": [
        {"source": "A", "relation": "binds", "target": "B"}
      ]
    }
  }
}
```

## 17) validate

- **作用：** 检查知识图谱的完整性：重复节点、悬空边、空的或不允许的关系、自环、重复边和孤立节点。

- **处理方法：** 错误：重复节点、指向未定义节点的边（给出了节点列表且 `allow_implicit_nodes` 为假时）、空关系、`allowed_relations` 之外的关系，以及 `allow_self_loops` 为假时的自环。警告：自环、重复边和孤立节点。最多列出 100 条错误和 100 条警告。只检查，不修复。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`allowed_relations`（列表；其他关系为错误）、`allow_self_loops`（`True`）、`allow_implicit_nodes`（`False`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `valid`（无错误）、`n_nodes`、`n_edges`、`errors` 与 `warnings`（含 `code`、`message`）。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import validate

graph = {"nodes": ["A", "B", "Z"], "edges": [["A", "binds", "B"], ["A", "binds", "B"], ["B", "binds", "B"],
                                            ["B", "binds", "C"]]}
result = validate(graph)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "validate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 121,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d8cd6b2142b90886f4328ea2c449ce225c76200d4c05bd0d359d628cab611f2d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "n_nodes": 4,
    "n_edges": 4,
    "errors": [
      {"code": "DANGLING_EDGE", "message": "edge 'B' -binds-> 'C' references undefined node 'C'"}
    ],
    "warnings": [
      {"code": "SELF_LOOP", "message": "self-loop on 'B' (binds)"},
      {"code": "DUPLICATE_EDGE", "message": "edge 'A' -binds-> 'B' appears 2 times"},
      {"code": "ISOLATED_NODE", "message": "1 node(s) have no edges, e.g. ['Z']"}
    ]
  }
}
```

## 18) version

- **作用：** 计算基于内容的版本号（规范化排序后图的 SHA-256），并给出节点数和边数。

- **处理方法：** 对节点、边及其属性的规范化（排序后）JSON 取 SHA-256，写作 `sha256:<hex>`；版本号与节点和边的顺序无关，任何内容变化都会改变它。给出 `previous_version` 时，还会报告内容自该版本号以来是否变化。纯 Python、本地运行；不使用图数据库、网络服务或可选库。

- **API：** `dnakit.knowledge_graph.version(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：图：映射 `{'nodes': [...], 'edges': [...], 'directed': bool}`、三元组列表 `[subject, relation, object]`、二元边或边映射（键 subject|source|head|from、predicate|relation|type、object|target|tail|to；其他键成为边的属性）、带分隔符的三元组文本，或 `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` 路径。通用选项：`directed`（默认 `True`，或取自输入映射）、`default_relation`（`related_to`）、`max_edges`（500000）、`drop_self_loops`（`False`）。`options`：`label`（原样返回）、`previous_version`（版本号；会增加 `changed_since_previous`）。

- **输出：** 标准封装（`status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings`）。`result` 含 `version_id`、`short_id`（前 12 位十六进制）、`directed`、`n_nodes`、`n_edges`、`relation_counts`、`definition`、原样返回的 `label`，给出 `previous_version` 时另有 `changed_since_previous`。

- **示例代码：**

```python
import json
from dnakit.knowledge_graph import version

result = version([["A", "binds", "B"], ["B", "binds", "C"]], options={"label": "release-1"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "knowledge_graph",
  "operation": "version",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 42,
    "sequence_bases": 14,
    "lengths": [7, 7],
    "sha256": "80f96dc5cc4579423529eb7d64469b6493b6959097d5721000531a1bc7f118c0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "version_id": "sha256:633b3250192822c51716b5bdaf9ba31af4ae8b886759d18da7760f649ec6c0d9",
    "short_id": "633b32501928",
    "directed": true,
    "n_nodes": 3,
    "n_edges": 2,
    "relation_counts": {"binds": 2},
    "definition": "SHA-256 of the canonical (sorted) JSON of nodes, edges and their properties",
    "label": "release-1"
  }
}
```
