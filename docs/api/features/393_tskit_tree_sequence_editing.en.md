# tskit tree-sequence editing and mutation tables

These functions in `dnakit.population` wrap the editing operations of [tskit](https://github.com/tskit-dev/tskit): table clean-up, interval and time cuts, trimming, simplification, subsetting, union and concatenation, edge splitting, haplotype extension and parsimony mutation mapping.

tskit is an optional dependency (`dnakit[bio-tskit]`) imported only when a function runs; a missing package raises an error with the install command. Every function takes the tree sequence as explicit table rows (the format of `build_tree_sequence()`) or as the path of a `.trees` file, and editing functions return the edited tables in the same row format. Intervals are half-open `[left, right)`. No network is used.

## Example data

All examples use one tree sequence of 4 samples and 100 bp: tree 1 on [0, 50) joins samples (0, 1) and (2, 3); tree 2 on [50, 100) joins (0, 2) and (1, 3). It has 4 sites with one mutation each.

**Example file:** `tables.json`

```text
{
 "sequence_length": 100,
 "nodes": [
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 0,
   "sample": true
  },
  {
   "time": 1
  },
  {
   "time": 2
  },
  {
   "time": 4
  },
  {
   "time": 1
  },
  {
   "time": 3
  },
  {
   "time": 5
  }
 ],
 "edges": [
  {
   "left": 0,
   "right": 50,
   "parent": 4,
   "child": 0
  },
  {
   "left": 0,
   "right": 50,
   "parent": 4,
   "child": 1
  },
  {
   "left": 0,
   "right": 50,
   "parent": 5,
   "child": 2
  },
  {
   "left": 0,
   "right": 50,
   "parent": 5,
   "child": 3
  },
  {
   "left": 0,
   "right": 50,
   "parent": 6,
   "child": 4
  },
  {
   "left": 0,
   "right": 50,
   "parent": 6,
   "child": 5
  },
  {
   "left": 50,
   "right": 100,
   "parent": 7,
   "child": 0
  },
  {
   "left": 50,
   "right": 100,
   "parent": 7,
   "child": 2
  },
  {
   "left": 50,
   "right": 100,
   "parent": 8,
   "child": 1
  },
  {
   "left": 50,
   "right": 100,
   "parent": 8,
   "child": 3
  },
  {
   "left": 50,
   "right": 100,
   "parent": 9,
   "child": 7
  },
  {
   "left": 50,
   "right": 100,
   "parent": 9,
   "child": 8
  }
 ],
 "sites": [
  {
   "position": 10,
   "ancestral_state": "A"
  },
  {
   "position": 30,
   "ancestral_state": "C"
  },
  {
   "position": 60,
   "ancestral_state": "G"
  },
  {
   "position": 80,
   "ancestral_state": "T"
  }
 ],
 "mutations": [
  {
   "site": 0,
   "node": 4,
   "derived_state": "T"
  },
  {
   "site": 1,
   "node": 5,
   "derived_state": "G"
  },
  {
   "site": 2,
   "node": 7,
   "derived_state": "A"
  },
  {
   "site": 3,
   "node": 8,
   "derived_state": "C"
  }
 ]
}
```

## 1) `dnakit.population.tskit_link_ancestors()`

- **Purpose:** Link the ancestral segments of chosen ancestor nodes to the samples that inherit them.
- **Method:** Calls tskit `TableCollection.link_ancestors`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_link_ancestors(tables, *, samples, ancestors)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `samples`, `ancestors` (node IDs).
- **Output:** `num_edges` and `edges` with `left`, `right`, `parent`, `child`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_link_ancestors

tables = json.load(open("tables.json"))
r = tskit_link_ancestors(tables, samples=[0, 1], ancestors=[4, 6])
print(r["num_edges"])
for e in r["edges"]:
    print(e)
```

- **Example output:**

```text
3
{'left': 0.0, 'right': 50.0, 'parent': 4, 'child': 0}
{'left': 0.0, 'right': 50.0, 'parent': 4, 'child': 1}
{'left': 0.0, 'right': 50.0, 'parent': 6, 'child': 4}
```

## 2) `dnakit.population.tskit_compute_mutation_parents()`

- **Purpose:** Compute the parent mutation of every mutation.
- **Method:** Calls tskit `TableCollection.compute_mutation_parents`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_compute_mutation_parents(tables)`.
- **Input:** `tables` (explicit rows or a `.trees` path).
- **Output:** `mutation_parents` (-1 means none) plus the edited tables.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_compute_mutation_parents

tables = json.load(open("tables.json"))
r = tskit_compute_mutation_parents(tables)
print(r["mutation_parents"])
```

- **Example output:**

```text
[-1, -1, -1, -1]
```

## 3) `dnakit.population.tskit_compute_mutation_times()`

- **Purpose:** Compute or impute mutation times from node times (`method="tables"` or `"min"`).
- **Method:** Calls tskit `TableCollection.compute_mutation_times` / `TreeSequence.impute_unknown_mutations_time`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_compute_mutation_times(tables, *, method='tables')`.
- **Input:** `tables` (explicit rows or a `.trees` path), `method`.
- **Output:** `mutation_times`; with `method="tables"` also the edited tables.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_compute_mutation_times

tables = json.load(open("tables.json"))
for method in ("tables", "min"):
    r = tskit_compute_mutation_times(tables, method=method)
    print(method, r["mutation_times"])
```

- **Example output:**

```text
tables [2.5, 3.0, 3.0, 4.0]
min [1.0, 2.0, 1.0, 3.0]
```

## 4) `dnakit.population.tskit_deduplicate_sites()`

- **Purpose:** Merge sites that lie at the same position.
- **Method:** Calls tskit `TableCollection.deduplicate_sites`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_deduplicate_sites(tables)`.
- **Input:** `tables` (explicit rows or a `.trees` path) (read as tables, because a tree sequence cannot hold duplicate sites).
- **Output:** `sites_before`, `sites_removed` and the edited tables.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_deduplicate_sites

tables = json.load(open("tables.json"))
dup = dict(tables, sites=tables["sites"] + [{"position": 10, "ancestral_state": "A"}], mutations=tables["mutations"] + [{"site": 4, "node": 1, "derived_state": "G"}])
r = tskit_deduplicate_sites(dup)
print(r["sites_before"], r["sites_removed"], r["result_summary"]["num_mutations"])
```

- **Example output:**

```text
5 1 5
```

## 5) `dnakit.population.tskit_delete_sites()`

- **Purpose:** Delete sites by ID together with their mutations.
- **Method:** Calls tskit `TreeSequence.delete_sites`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_delete_sites(tables, *, site_ids)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `site_ids`.
- **Output:** `deleted` and the edited tables.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_delete_sites

tables = json.load(open("tables.json"))
r = tskit_delete_sites(tables, site_ids=[0, 3])
print(r["deleted"], r["result_summary"]["num_sites"], [s["position"] for s in r["tables"]["sites"]])
```

- **Example output:**

```text
[0, 3] 2 [30.0, 60.0]
```

## 6) `dnakit.population.tskit_select_intervals()`

- **Purpose:** Delete or keep genomic intervals (`action="delete"` or `"keep"`).
- **Method:** Calls tskit `TreeSequence.delete_intervals` / `keep_intervals`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_select_intervals(tables, *, intervals, action='delete', simplify=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `intervals` as [left, right] pairs, `action`, `simplify`.
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_select_intervals

tables = json.load(open("tables.json"))
for action in ("delete", "keep"):
    r = tskit_select_intervals(tables, intervals=[[0, 50]], action=action)
    print(action, r["result_summary"]["num_trees"], r["result_summary"]["num_sites"])
```

- **Example output:**

```text
delete 2 2
keep 2 2
```

## 7) `dnakit.population.tskit_remove_older_genealogy()`

- **Purpose:** Remove the genealogy older than a time, with or without new roots.
- **Method:** Calls tskit `TreeSequence.decapitate` / `TableCollection.delete_older`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_remove_older_genealogy(tables, *, time, add_roots=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `time`, `add_roots`.
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_remove_older_genealogy

tables = json.load(open("tables.json"))
r = tskit_remove_older_genealogy(tables, time=3.5)
print(r["result_summary"]["num_edges"], r["result_summary"]["num_nodes"])
r = tskit_remove_older_genealogy(tables, time=3.5, add_roots=False)
print(r["result_summary"]["num_edges"], r["result_summary"]["num_nodes"])
```

- **Example output:**

```text
12 14
8 10
```

## 8) `dnakit.population.tskit_trim_genome()`

- **Purpose:** Trim empty genome flanks on the left, right or both sides.
- **Method:** Calls tskit `TreeSequence.ltrim` / `rtrim` / `trim`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_trim_genome(tables, *, side='both')`.
- **Input:** `tables` (explicit rows or a `.trees` path), `side`.
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_trim_genome

tables = json.load(open("tables.json"))
from dnakit.population.tskit_p31_functions import tskit_shift_coordinates

shifted = tskit_shift_coordinates(tables, value=10)["tables"]
print(shifted["sequence_length"], shifted["edges"][0]["left"])
r = tskit_trim_genome(shifted, side="left")
print(r["tables"]["sequence_length"], r["tables"]["edges"][0]["left"])
```

- **Example output:**

```text
110.0 10.0
100.0 0.0
```

## 9) `dnakit.population.tskit_shift_coordinates()`

- **Purpose:** Shift all genomic coordinates by a constant.
- **Method:** Calls tskit `TreeSequence.shift`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_shift_coordinates(tables, *, value, sequence_length=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `value`, `sequence_length`.
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_shift_coordinates

tables = json.load(open("tables.json"))
r = tskit_shift_coordinates(tables, value=10)
print(r["tables"]["sequence_length"], r["tables"]["sites"][0]["position"])
```

- **Example output:**

```text
110.0 20.0
```

## 10) `dnakit.population.tskit_simplify()`

- **Purpose:** Simplify the genealogy to the history of chosen samples.
- **Method:** Calls tskit `TreeSequence.simplify`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_simplify(tables, *, samples=None, reduce_to_site_topology=False, keep_unary=False, keep_input_roots=False, filter_sites=True, filter_nodes=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `samples` and the tskit simplify flags.
- **Output:** The edited tables, `node_map` and `options`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_simplify

tables = json.load(open("tables.json"))
r = tskit_simplify(tables, samples=[0, 1])
print(r["result_summary"]["num_nodes"], r["node_map"])
```

- **Example output:**

```text
4 [0, 1, -1, -1, 2, -1, -1, -1, -1, 3]
```

## 11) `dnakit.population.tskit_subset()`

- **Purpose:** Keep only the listed nodes and what lies between them.
- **Method:** Calls tskit `TreeSequence.subset`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_subset(tables, *, nodes, reorder_populations=True, remove_unreferenced=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `nodes`.
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_subset

tables = json.load(open("tables.json"))
r = tskit_subset(tables, nodes=[0, 1, 4])
print(r["result_summary"]["num_nodes"], r["result_summary"]["num_edges"])
```

- **Example output:**

```text
3 2
```

## 12) `dnakit.population.tskit_union()`

- **Purpose:** Add the non-shared part of a second tree sequence, matched node by node.
- **Method:** Calls tskit `TreeSequence.union`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_union(tables, other, *, node_mapping, all_edges=False, all_mutations=False, check_shared_equality=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `other`, `node_mapping` (node of `tables` per node of `other`, -1 for new).
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_union

tables = json.load(open("tables.json"))
r = tskit_union(tables, tables, node_mapping=list(range(10)))
print(r["result_summary"]["num_nodes"], r["result_summary"]["num_edges"])
```

- **Example output:**

```text
10 12
```

## 13) `dnakit.population.tskit_concatenate()`

- **Purpose:** Join tree sequences end to end along the genome.
- **Method:** Calls tskit `TreeSequence.concatenate`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_concatenate(tables, others)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `others` (a list of tree sequences).
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_concatenate

tables = json.load(open("tables.json"))
r = tskit_concatenate(tables, [tables])
print(r["result_summary"]["sequence_length"], r["result_summary"]["num_trees"], r["result_summary"]["num_sites"])
```

- **Example output:**

```text
200.0 4 8
```

## 14) `dnakit.population.tskit_extend_haplotypes()`

- **Purpose:** Extend ancestral haplotype segments through unary regions.
- **Method:** Calls tskit `TreeSequence.extend_haplotypes`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_extend_haplotypes(tables, *, max_iter=10)`.
- **Input:** `tables` (explicit rows or a `.trees` path) with known mutation times, `max_iter`.
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_extend_haplotypes

tables = json.load(open("tables.json"))
from dnakit.population.tskit_p31_functions import tskit_compute_mutation_times

timed = tskit_compute_mutation_times(tables)["tables"]
r = tskit_extend_haplotypes(timed)
print(r["result_summary"]["num_edges"], r["result_summary"]["num_trees"])
```

- **Example output:**

```text
12 2
```

## 15) `dnakit.population.tskit_split_edges()`

- **Purpose:** Insert a unary node on every edge that spans a given time.
- **Method:** Calls tskit `TreeSequence.split_edges`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_split_edges(tables, *, time)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `time`.
- **Output:** `result_summary` of the edited tree sequence and its `tables` in the input row format, plus `provenance`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_split_edges

tables = json.load(open("tables.json"))
r = tskit_split_edges(tables, time=1.5)
print(r["result_summary"]["num_nodes"], r["result_summary"]["num_edges"])
```

- **Example output:**

```text
16 18
```

## 16) `dnakit.population.tskit_map_mutations()`

- **Purpose:** Place mutations on one tree that explain observed genotypes by parsimony.
- **Method:** Calls tskit `Tree.map_mutations`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_map_mutations(tables, *, genotypes, tree_index=0, alleles=None, ancestral_state=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `genotypes` (one allele index per sample), `tree_index`, `alleles`, `ancestral_state`.
- **Output:** `ancestral_state`, `alleles` and the `mutations` (node, derived_state, parent, time).
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_map_mutations

tables = json.load(open("tables.json"))
r = tskit_map_mutations(tables, genotypes=[1, 1, 0, 0], tree_index=0)
print(r["ancestral_state"], r["mutations"])
```

- **Example output:**

```text
0 [{'node': 4, 'derived_state': '1', 'parent': -1, 'time': None}]
```
