# Knowledge graph

This category contains **18** documented functions. 17 of them operate on graphs held as plain triples/edge lists in pure Python (no graph database or web service): all 17 compute their result locally and none of them returns `not_implemented`. The remaining one, `pseudonymize_genomic_metadata_fields`, is unrelated to graphs: it pseudonymises sensitive fields of a metadata mapping (key-name and pattern based, not de-identification) and is described in its own section. Graph inputs can be a mapping, a list of triples, delimited text or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path, so all functions can be called through Python and through JSON/MCP. Note the default of `apply_access_policy`: items without a classification label are removed unless `unclassified='allow'` is given. Examples below shorten long arrays and text with `...`.

## 1) add_graph_relationship

- **Purpose:** Add one or more relationships (edges) to a knowledge graph and return the new graph; the input is not modified.

- **Method:** Loads the graph (or starts from an empty one when the input is `None`), adds the given edge(s) with their properties and creates missing nodes unless `create_missing_nodes` is false. An edge is identified by (source, relation, target); a repeated edge is handled by `on_duplicate`. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.add_graph_relationship(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`); or `None`. `options`: `subject`, `predicate`, `object` and optional `properties` for one edge, or `relationships` (list of triples / edge mappings); `create_missing_nodes` (`True`; if `False` unknown nodes raise), `on_duplicate` (`ignore` default, `error`, `update`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `added` (new edges), `already_present` (edges that existed) and `graph` (`directed`, `nodes`, `edges`).

- **Example code:**

```python
import json
from dnakit.knowledge_graph import add_graph_relationship

graph = [["BRCA1", "interacts_with", "TP53"]]
result = add_graph_relationship(graph, options={"subject": "TP53", "predicate": "regulates",
                                                "object": "CDKN1A", "properties": {"source": "demo"}})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Aggregate a knowledge graph into summary statistics: relation counts, degree statistics, density, components, node types and top hubs.

- **Method:** Degree = in-edges + out-edges of the graph's edge set (a self-loop counts twice); density = edges / possible edges excluding self-loops; components are weakly connected components; `top_hubs` are the `top_k` nodes of highest degree; node types are counted from the `type_property` of each node (`unknown` when absent). Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.aggregate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `top_k` (10), `type_property` (`type`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `degree` (`total`, `in`, `out`, each with min/max/mean/median), `density`, `n_isolated_nodes`, `n_self_loops`, `n_weakly_connected_components`, `largest_component_size`, `node_type_counts`, `top_hubs` and a `degree_note`.

- **Example code:**

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

- **Example output:**

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

- **Purpose:** Apply a clearance-level access policy to a knowledge graph: remove nodes and edges classified above the clearance and redact properties.

- **Method:** Each node and edge carries a label in the property `classification_property`; the `levels` list orders the labels from low to high. An item is kept when its level is at or below `clearance`; edges whose endpoint node was removed are removed too, then the properties in `redact_properties` are dropped from the kept items. **Defaults to know about:** `clearance` defaults to the lowest level (`public`), and `unclassified` defaults to `deny`, i.e. every node or edge without a classification label is removed. A graph with no classification labels, called with default options, therefore comes back empty. Pass `unclassified='allow'` to keep unclassified items. A label that is not in `levels` raises `UNKNOWN_CLASSIFICATION`. This filters data; it is not an authentication system. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.apply_access_policy(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `clearance` (default the lowest level), `levels` (low to high; default `['public', 'internal', 'restricted', 'secret']`), `classification_property` (`classification`, read from node and edge properties), `unclassified` (`deny` default: items without a classification are removed; or `allow`), `deny_nodes` (ids removed regardless), `redact_properties` (property names dropped from kept items).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts` of the filtered graph, `clearance`, `levels`, `n_nodes_removed`, `removed_nodes` (at most 1000), `n_edges_removed`, a `policy` sentence stating whether unclassified items were removed or kept, and `graph`.

- **Example code:**

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

- **Example output:**

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

- **Purpose:** Build a data-lineage graph (output derived_from input) from provenance records and report roots, leaves and a topological order.

- **Method:** Every record contributes edges from its output to each of its inputs with relation `derived_from`, the step name stored on the edge. Roots are the nodes with no inputs, leaves the nodes nothing is derived from; the topological order lists inputs before outputs. A cycle raises an error. With `trace`, all upstream ancestors of that node are returned with their distance (breadth-first). Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.build_data_lineage(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of records `{'output': id (or 'target' / 'dataset' / 'id'), 'inputs': [ids] (or 'sources' / 'derived_from'), 'step': name (or 'tool' / 'process')}`. `options`: `trace` (node id).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `roots`, `leaves`, `topological_order`, `edge_meaning`, `graph` and, with `trace`, `upstream_of` (`node`, `ancestors` with `id` and `distance`).

- **Example code:**

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

- **Example output:**

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

- **Purpose:** Build a knowledge graph (unique nodes, merged duplicate edges) from triples or an edge list.

- **Method:** Nodes mentioned only by edges are created; repeated (source, relation, target) edges are merged and the number of repetitions is stored in the edge property `count`. Self-loops are kept unless `drop_self_loops` is true. For undirected graphs (`directed=False`) edge endpoints are treated as unordered. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.build_knowledge_graph(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `duplicate_nodes`, `duplicate_edges`, `self_loops`, `dropped_self_loops` and `graph` (`directed`, `nodes`, `edges`).

- **Example code:**

```python
import json
from dnakit.knowledge_graph import build_knowledge_graph

triples = [["A", "binds", "B"], ["A", "binds", "B"], ["B", "binds", "B"], ["C", "binds", "A"]]
result = build_knowledge_graph(triples)
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Compare two versions of a knowledge graph: added, removed and changed nodes and edges.

- **Method:** Nodes are identified by id and edges by (source, relation, target); a node or edge is "changed" when its properties differ. The two content-based version ids (see `version`) are reported and `identical` is true when they are equal. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.compare_database_versions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'old': graph, 'new': graph}`, a pair `[old, new]` of graph mappings or files, or one graph plus `options['other']` as the new version. Graph forms as in `build_knowledge_graph`. `options`: `max_items` (1000, bounds each list).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `identical`, `old_version_id`, `new_version_id`, `counts` (nodes/edges added, removed, changed), `added_nodes`, `removed_nodes`, `changed_nodes` (`id`, `old`, `new` properties), `added_edges`, `removed_edges`, `changed_edges` and `truncated`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import compare_database_versions

old = {"nodes": [{"id": "A", "score": 1}, "B"], "edges": [["A", "binds", "B"]]}
new = {"nodes": [{"id": "A", "score": 2}, "B", "C"], "edges": [["A", "binds", "B"], ["B", "binds", "C"]]}
result = compare_database_versions({"old": old, "new": new})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Create a new knowledge graph, optionally seeded with nodes and edges.

- **Method:** Starts from an empty graph (input `None`) or from any graph input accepted by `build_knowledge_graph`, then adds the `nodes` and `edges` options and stores `name` and `description` as graph metadata. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.create(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `None` for an empty graph, or any graph input accepted by `build_knowledge_graph`. `options`: `name`, `description`, `directed` (`True`), `nodes` (list of ids or `{'id': ..., properties}`), `edges` (list of triples or edge mappings), `default_relation`.

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `duplicate_nodes`, `duplicate_edges`, `self_loops`, `dropped_self_loops` and `graph` (with `name` and `description` when given).

- **Example code:**

```python
import json
from dnakit.knowledge_graph import create

result = create(None, options={"name": "demo", "nodes": ["A", {"id": "B", "type": "gene"}],
                               "edges": [["A", "binds", "B"]]})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Pseudonymise the sensitive fields of a metadata mapping (or list of mappings): replace, drop or keyed-hash the values of identifying keys at any nesting depth, and replace e-mail addresses, phone numbers and SSN-like patterns inside text. This is key-name and pattern based pseudonymisation, NOT a guarantee of de-identification. It replaces an earlier function that hashed only six top-level keys with an unsalted SHA-256; that function was removed because its name promised more than it did.

- **Method:** Implemented once in `src/dnakit/_metadata_pseudonymize.py` and shared with `dnakit.workflow_security.pseudonymize_metadata_fields`. The input is never modified; a deep copy is transformed. A key is sensitive when its name, lower-cased with everything except letters and digits removed (so `Patient_ID`, `patient id` and `PATIENTID` are the same), is in the default list: `name`, `first/last/full/given/family/middle name`, `surname`, `patient`, `patient id`, `patient name`, `mrn`, `medical record number`, `ssn`, `social security number`, `email`, `email address`, `phone`, `phone number`, `telephone`, `mobile`, `fax`, `address`, `home/street address`, `street`, `city`, `zip`, `zip code`, `postal code`, `postcode`, `dob`, `date of birth`, `birth date`, `birthday`, `national id`, `passport`, `passport number`, `insurance`, `insurance id/number`, `account number`, `ip address`, `device id`, `photo`, `photograph`, `biometric(s)`, `license number`, `driver license`. Matching is exact on the normalised name (`sample_name` is not matched); `options['keys']` adds names and `options['keep']` removes names (keep wins). The value of a sensitive key is replaced as a whole, also when it is a mapping or list. `mode`: `redact` (default) writes `"[REDACTED]"`, needing no secret; `drop` removes the key; `hash` writes the first `digest_chars` (default 16, 8 to 64) hexadecimal characters of HMAC-SHA256 keyed with `options['salt']` (required, a string of at least 8 characters, otherwise `ConfigurationError`) over the UTF-8 value (strings as they are, other values as compact sorted-key JSON). Hash tokens are deterministic for one salt, so equal values stay linkable and low-entropy values can be guessed by anyone who knows the salt; keep the salt secret. With `scan_text` (default `True`) every string value anywhere, including values under kept keys and in lists, is scanned for e-mail addresses, SSN-like numbers (`123-45-6789`), international phone numbers (`+44 20 7946 0958`) and North-American style 10-digit numbers; matches become `[REDACTED:email]`, `[REDACTED:ssn]`, `[REDACTED:phone]`. The phone pattern also matches bare 10-digit identifiers, so it can over-redact. What it does NOT do: dates are not shifted or generalised, ages over 89 and small geographic units are not generalised, names inside free text and other quasi-identifiers are not detected, dictionary keys that themselves contain identifiers are not rewritten, and nothing is checked against HIPAA or GDPR identifier lists. Review the output before sharing it. Datetimes become ISO strings; other non-JSON types (sets, bytes, objects) and nesting deeper than 64 levels or cycles raise `InputFormatError`.

- **API:** `dnakit.knowledge_graph.pseudonymize_genomic_metadata_fields(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a mapping, a list of mappings, or JSON text of one of them; anything else raises `InputFormatError` (codes `INVALID_METADATA_INPUT`, `INVALID_METADATA_JSON`, `UNSUPPORTED_METADATA_VALUE`, `METADATA_TOO_DEEP`, each with a `hint`). `options`: `mode` (`"redact"` | `"drop"` | `"hash"`, default `"redact"`), `salt` (string of at least 8 characters, required for `hash`), `digest_chars` (8 to 64, default 16, `hash` only), `keys` (extra sensitive key names, a string or list), `keep` (key names never treated as sensitive), `scan_text` (default `True`), `backend` (label only). Invalid options raise `ConfigurationError` (`PSEUDONYMIZE_SALT_REQUIRED`, `INVALID_PSEUDONYMIZE_MODE`, `INVALID_DIGEST_CHARS`, `INVALID_PSEUDONYMIZE_KEYS`, `INVALID_PSEUDONYMIZE_OPTIONS`).

- **Output:** The standard envelope: `status` (`"ok"`), `module` (`"knowledge_graph"`), `operation` (`"pseudonymize_genomic_metadata_fields"`), `provenance` (`implementation` `"native"`, `backend`), `input_summary`, `warnings` and `result`. `input_summary` deliberately holds only `record_count` and `n_keys` (number of mapping keys visited): no length and no SHA-256 of the unprotected input. `result` has `metadata` (the transformed mapping or list; JSON text input is returned parsed), `mode`, `n_redacted` (= `n_key_matches` + `n_text_matches`), `n_key_matches`, `n_text_matches`, `text_matches_by_pattern`, `redacted_paths` (JSON-path-like key paths such as `$.visits[0].dob`, never values; at most 10000, `redacted_paths_truncated` tells if more existed) and `warnings`. The warnings always say that the method is key-name and pattern based, is NOT a guarantee of de-identification, does not shift or generalise dates, ages over 89 or small geographic units, that free text may still identify a person, and that the output must be reviewed before sharing; the same list is in `result.warnings` and the envelope `warnings`. Neither the salt nor any original value appears anywhere in the output.

- **Example code:**

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

- **Example output:**

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

- **Purpose:** Serialise a knowledge graph as JSON, JSONL, TSV (edge triples only) or GraphML text, optionally writing it to a file.

- **Method:** GraphML stores property values as strings; TSV holds the edge triples only, so node and edge properties and isolated nodes are dropped (the result has a `note`). The file is written atomically and an existing file is not replaced unless `overwrite` is true. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `format` (`json`, `jsonl`, `tsv`, `graphml`; default from the `output` suffix, else `json`), `output` (path; without it the text is returned under `content`), `overwrite` (`False`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `format`, `bytes`, `content` (when no `output` is given), `output_path` and, for TSV, a `note`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import export

result = export([["A", "binds", "B"], ["B", "regulates", "C"]], options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Link entities: create edges between explicit node pairs, or between all nodes that share the same value of a node property.

- **Method:** With `pairs`, an edge with relation `relation` is added for each [a, b] (missing nodes are created). With `by`, nodes whose value of that property is equal (ignoring case) are linked pairwise, for groups of at most `max_group_size` nodes. Edges that already exist are not duplicated. Exactly one of `pairs` or `by` must be given. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.link(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `relation` (`linked_to`); exactly one of `pairs` (list of `[a, b]`) or `by` (node property name); `max_group_size` (100).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `n_created`, `created` (the new edges) and `graph`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import link

graph = {"nodes": [{"id": "g1", "pathway": "p53"}, {"id": "g2", "pathway": "P53"}, {"id": "g3", "pathway": "wnt"}],
         "edges": []}
result = link(graph, options={"by": "pathway", "relation": "same_pathway"})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Map (rename) node identifiers and relation names through a lookup table, merging nodes that map to the same identifier.

- **Method:** Node ids are replaced through `mapping` and relation names through `relation_mapping`. Nodes that map to the same new id are merged (the merged node keeps the properties of the first node) and the edges that become identical are merged and counted in `count`. Nodes absent from `mapping` are kept, dropped or rejected according to `unmapped`. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `mapping` (dict old id -> new id, rows `{'from', 'to'}` / `[from, to]`, or a two-column TSV/CSV path), `relation_mapping` (same forms, for relation names), `unmapped` (`keep` default, `drop`, `error`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `n_ids_renamed`, `n_nodes_merged`, `n_nodes_dropped` and `graph`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import map as map_graph

graph = [["ENSG1", "interacts_with", "ENSG2"], ["ENSG3", "interacts_with", "ENSG2"]]
result = map_graph(graph, options={"mapping": {"ENSG1": "BRCA1", "ENSG2": "TP53", "ENSG3": "BRCA1"}})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Query the neighbourhood of a node by breadth-first search, optionally with the shortest path to a target and the induced subgraph.

- **Method:** Breadth-first search from `node` up to `depth` steps, following edges in the given `direction` (ignored for undirected graphs) and only the listed `relations` if given; each neighbour is reported with its distance, one shortest path and the edge it was reached through. With `target`, the shortest path (fewest edges, under the same direction and relation rules, at most `max_path_length`) is returned. A unknown `node` raises an error. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.query_graph_neighbors(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `node` (required), `depth` (1, at most 20), `direction` (`both` default, `out`, `in`), `relations` (restrict to these relation names), `limit` (1000), `target`, `max_path_length` (50), `include_subgraph` (`False`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `node`, `depth`, `direction`, `n_neighbors`, `neighbors` (`id`, `distance`, `path`, `via_edge`), `truncated`, and with `target` a `shortest_path` (`found`, `length`, `nodes`, `edges`), and with `include_subgraph` the induced `subgraph`.

- **Example code:**

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

- **Example output:**

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

- **Purpose:** Read a knowledge graph from a file or text (JSON, JSONL, TSV/CSV triples, GraphML) and return it with its summary.

- **Method:** Parses the input as described for the graph forms and builds the graph like `build_knowledge_graph`. GraphML with DOCTYPE/ENTITY declarations is rejected, and files larger than 100 MB are refused. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `duplicate_nodes`, `duplicate_edges`, `self_loops`, `dropped_self_loops` and `graph`.

- **Example code:**

```python
import json
from pathlib import Path
from dnakit.knowledge_graph import read

Path("edges.tsv").write_text("BRCA1\tinteracts_with\tTP53\nTP53\tregulates\tCDKN1A\n")
result = read("edges.tsv")
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Search the nodes of a knowledge graph by text in their id and string or number properties.

- **Method:** Compares the query with the node id and with every string or number property (restricted to `fields` if given) using `match` (substring, exact or prefix; case-insensitive unless `case_sensitive`). Nodes can be restricted to one `node_type`. Matches keep the graph's node order; the fields that matched are reported. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.search(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `query` (required), `match` (`contains` default, `exact`, `prefix`), `case_sensitive` (`False`), `fields` (property names to search; `id` included), `node_type` with `type_property` (`type`), `limit` (100).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `query`, `match`, `n_matches`, `matches` (node `id`, `matched_fields` and the node's properties) and `truncated`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import search

graph = {"nodes": [{"id": "BRCA1", "type": "gene", "description": "DNA repair"},
                   {"id": "TP53", "type": "gene", "description": "tumour suppressor"},
                   {"id": "repair_pathway", "type": "pathway"}], "edges": []}
result = search(graph, options={"query": "repair"})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Create a snapshot of a knowledge graph (canonical content plus version id), optionally saved as a JSON file.

- **Method:** The snapshot is the canonical (sorted) graph with its content-based `version_id` (as in `version`), the optional `label` and, only if given, the `timestamp`, so repeated snapshots of the same graph are byte-identical. The file is written atomically and an existing file is not replaced unless `overwrite` is true. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.snapshot_database_version(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `label`, `timestamp` (stored only if given), `output` (path; otherwise the snapshot is returned), `overwrite` (`False`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `version_id`, `directed`, `n_nodes`, `n_edges`, `relation_counts`, `snapshot` (`format`, `format_version`, `version_id`, `label`, `graph`), `bytes`, `content` (JSON text, when no `output` is given) and `output_path`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import snapshot_database_version

result = snapshot_database_version([["A", "binds", "B"]], options={"label": "v1"})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Apply a change set to a knowledge graph: remove edges or nodes (with their edges) and set node or edge properties.

- **Method:** Returns a new graph with the requested changes applied; the input is not modified. Removing a node removes all its edges. Unknown ids raise an error unless `ignore_missing` is true. At least one change option is required. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.update(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options` (at least one required): `remove_edges` (list of `[source, relation, target]`), `remove_nodes` (ids), `set_node_properties` (`{id: {property: value}}`), `set_edge_properties` (list of `{source, relation, target, properties}`), `ignore_missing` (`False`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `directed`, `n_nodes`, `n_edges` and `relation_counts`, `edges_removed`, `nodes_removed`, `node_properties_set`, `edge_properties_set` and `graph`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import update

graph = [["A", "binds", "B"], ["B", "binds", "C"]]
result = update(graph, options={"remove_nodes": ["C"], "set_node_properties": {"A": {"type": "gene"}}})
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Validate the integrity of a knowledge graph: duplicate nodes, dangling edges, empty or disallowed relations, self-loops, duplicate edges and isolated nodes.

- **Method:** Errors: duplicate nodes, edges to undefined nodes (when a node list is given and `allow_implicit_nodes` is false), empty relations, relations outside `allowed_relations`, and self-loops when `allow_self_loops` is false. Warnings: self-loops, duplicate edges and isolated nodes. At most 100 errors and 100 warnings are listed. The graph is only inspected, never repaired. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `allowed_relations` (list; others are errors), `allow_self_loops` (`True`), `allow_implicit_nodes` (`False`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `valid` (no errors), `n_nodes`, `n_edges`, `errors` and `warnings` (items with `code`, `message`).

- **Example code:**

```python
import json
from dnakit.knowledge_graph import validate

graph = {"nodes": ["A", "B", "Z"], "edges": [["A", "binds", "B"], ["A", "binds", "B"], ["B", "binds", "B"],
                                            ["B", "binds", "C"]]}
result = validate(graph)
print(json.dumps(result, indent=2))
```

- **Example output:**

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

- **Purpose:** Compute a content-based version identifier (SHA-256 of the canonical sorted graph) with node and edge counts.

- **Method:** SHA-256 of the canonical (sorted) JSON of nodes, edges and their properties, written as `sha256:<hex>`; the id does not depend on node or edge order, and any change of content changes it. With `previous_version` the call also reports whether the content changed since that id. Pure-Python, local; no graph database, web service or optional library is used.

- **API:** `dnakit.knowledge_graph.version(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a graph: a mapping `{'nodes': [...], 'edges': [...], 'directed': bool}`, a list of triples `[subject, relation, object]`, 2-item edges or edge mappings (keys subject|source|head|from, predicate|relation|type, object|target|tail|to; other keys become edge properties), delimited triple text, or a `.json`/`.jsonl`/`.tsv`/`.csv`/`.graphml` path. Common options: `directed` (`True`, or taken from the input mapping), `default_relation` (`related_to`), `max_edges` (500000), `drop_self_loops` (`False`). `options`: `label` (echoed), `previous_version` (a version id; adds `changed_since_previous`).

- **Output:** The standard envelope (`status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings`). `result` has `version_id`, `short_id` (first 12 hex digits), `directed`, `n_nodes`, `n_edges`, `relation_counts`, `definition`, the echoed `label` and, with `previous_version`, `changed_since_previous`.

- **Example code:**

```python
import json
from dnakit.knowledge_graph import version

result = version([["A", "binds", "B"], ["B", "binds", "C"]], options={"label": "release-1"})
print(json.dumps(result, indent=2))
```

- **Example output:**

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
