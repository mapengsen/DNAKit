# tskit listings, genotypes, sequences, IBD and recombination maps

These functions in `dnakit.population` wrap the listing and export operations of [tskit](https://github.com/tskit-dev/tskit): mutations, sites, genotype matrices, variants, alignments, haplotypes, breakpoints, VCF individual mapping, IBD segments and recombination-map look-ups.

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

**Example file:** `map.txt`

```text
Chromosome	Position(bp)	Rate(cM/Mb)	Map(cM)
chr1	0	1.0	0.0
chr1	1000000	2.0	1.0
chr1	2000000	0.0	3.0
```

## 1) `dnakit.population.tskit_list_mutations()`

- **Purpose:** List the mutations of the tree sequence or of one tree.
- **Method:** Calls tskit `TreeSequence.mutations` / `Tree.mutations`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_list_mutations(tables, *, tree_index=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `tree_index`.
- **Output:** `num_mutations` and `mutations` rows.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_list_mutations

tables = json.load(open("tables.json"))
r = tskit_list_mutations(tables, tree_index=1)
print(r["num_mutations"])
for m in r["mutations"]:
    print(m)
```

- **Example output:**

```text
2
{'id': 2, 'site': 2, 'node': 7, 'derived_state': 'A', 'parent': -1, 'time': None, 'position': 60.0}
{'id': 3, 'site': 3, 'node': 8, 'derived_state': 'C', 'parent': -1, 'time': None, 'position': 80.0}
```

## 2) `dnakit.population.tskit_list_sites()`

- **Purpose:** List the variant sites of the tree sequence or of one tree.
- **Method:** Calls tskit `TreeSequence.sites` / `Tree.sites`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_list_sites(tables, *, tree_index=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `tree_index`.
- **Output:** `num_sites` and `sites` rows.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_list_sites

tables = json.load(open("tables.json"))
r = tskit_list_sites(tables)
for s in r["sites"]:
    print(s)
```

- **Example output:**

```text
{'id': 0, 'position': 10.0, 'ancestral_state': 'A', 'num_mutations': 1}
{'id': 1, 'position': 30.0, 'ancestral_state': 'C', 'num_mutations': 1}
{'id': 2, 'position': 60.0, 'ancestral_state': 'G', 'num_mutations': 1}
{'id': 3, 'position': 80.0, 'ancestral_state': 'T', 'num_mutations': 1}
```

## 3) `dnakit.population.tskit_genotype_matrix()`

- **Purpose:** Return the genotype matrix of allele indexes (sites by samples).
- **Method:** Calls tskit `TreeSequence.genotype_matrix`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_genotype_matrix(tables, *, samples=None, isolated_as_missing=True, alleles=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `samples`, `isolated_as_missing`, `alleles`.
- **Output:** `matrix` (-1 is missing), `positions`, `num_sites`, `num_samples`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genotype_matrix

tables = json.load(open("tables.json"))
r = tskit_genotype_matrix(tables)
print(r["positions"])
for row in r["matrix"]:
    print(row)
```

- **Example output:**

```text
[10.0, 30.0, 60.0, 80.0]
[1, 1, 0, 0]
[0, 0, 1, 1]
[1, 0, 1, 0]
[0, 1, 0, 1]
```

## 4) `dnakit.population.tskit_variants()`

- **Purpose:** Per-site genotype calls with allele states and allele counts.
- **Method:** Calls tskit `TreeSequence.variants`, `Variant.states`, `Variant.counts`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_variants(tables, *, samples=None, isolated_as_missing=True, alleles=None, left=None, right=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `samples`, `isolated_as_missing`, `alleles`, `left`, `right`.
- **Output:** `variants` rows with `site`, `position`, `alleles`, `genotypes`, `states`, `counts`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_variants

tables = json.load(open("tables.json"))
r = tskit_variants(tables)
for v in r["variants"][:2]:
    print(v["position"], v["alleles"], v["genotypes"], v["states"], v["counts"])
```

- **Example output:**

```text
10.0 ['A', 'T'] [1, 1, 0, 0] [np.str_('T'), np.str_('T'), np.str_('A'), np.str_('A')] {'A': 2, 'T': 2}
30.0 ['C', 'G'] [0, 0, 1, 1] [np.str_('C'), np.str_('C'), np.str_('G'), np.str_('G')] {'C': 2, 'G': 2}
```

## 5) `dnakit.population.tskit_alignments()`

- **Purpose:** Full-length sample sequence alignments over a reference.
- **Method:** Calls tskit `TreeSequence.alignments`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_alignments(tables, *, reference_sequence, missing_data_character='N', isolated_as_missing=True, samples=None, left=None, right=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `reference_sequence` (length = sequence length), `missing_data_character`, `samples`, `left`, `right`.
- **Output:** `alignments` rows with `node` and `sequence`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_alignments

tables = json.load(open("tables.json"))
r = tskit_alignments(tables, reference_sequence="A" * 100)
for a in r["alignments"]:
    print(a["node"], a["sequence"][8:12], a["sequence"][58:62])
```

- **Example output:**

```text
0 AATA AAAA
1 AATA AAGA
2 AAAA AAAA
3 AAAA AAGA
```

## 6) `dnakit.population.tskit_haplotypes()`

- **Purpose:** Haplotype string per sample over the variable sites only.
- **Method:** Calls tskit `TreeSequence.haplotypes`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_haplotypes(tables, *, missing_data_character='-', isolated_as_missing=True, samples=None, left=None, right=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `missing_data_character`, `isolated_as_missing`, `samples`, `left`, `right`.
- **Output:** `haplotypes` rows with `node` and `haplotype`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_haplotypes

tables = json.load(open("tables.json"))
r = tskit_haplotypes(tables)
for h in r["haplotypes"]:
    print(h["node"], h["haplotype"])
```

- **Example output:**

```text
0 TCAT
1 TCGC
2 AGAT
3 AGGC
```

## 7) `dnakit.population.tskit_breakpoints()`

- **Purpose:** Genomic breakpoints between genealogies.
- **Method:** Calls tskit `TreeSequence.breakpoints`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_breakpoints(tables)`.
- **Input:** `tables` (explicit rows or a `.trees` path).
- **Output:** `breakpoints` (starting at 0, ending at the sequence length) and `num_trees`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_breakpoints

tables = json.load(open("tables.json"))
r = tskit_breakpoints(tables)
print(r["breakpoints"], r["num_trees"])
```

- **Example output:**

```text
[0.0, 50.0, 100.0] 2
```

## 8) `dnakit.population.tskit_map_to_vcf_model()`

- **Purpose:** Map sample nodes to VCF individuals.
- **Method:** Calls tskit `TreeSequence.map_to_vcf_model`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_map_to_vcf_model(tables, *, ploidy=None, contig_id='1')`.
- **Input:** `tables` (explicit rows or a `.trees` path), `ploidy`, `contig_id`.
- **Output:** `individual_names`, `individual_nodes`, `contig_id`, `contig_length`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_map_to_vcf_model

tables = json.load(open("tables.json"))
r = tskit_map_to_vcf_model(tables, ploidy=2)
print(r["individual_names"], r["individual_nodes"], r["contig_length"])
```

- **Example output:**

```text
['tsk_0', 'tsk_1'] [[0, 1], [2, 3]] 100
```

## 9) `dnakit.population.tskit_sample_nodes_by_ploidy()`

- **Purpose:** Group the sample nodes into rows of a given ploidy.
- **Method:** Calls tskit `TreeSequence.sample_nodes_by_ploidy`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_sample_nodes_by_ploidy(tables, *, ploidy)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `ploidy`.
- **Output:** `groups`, one list of node IDs per individual.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_sample_nodes_by_ploidy

tables = json.load(open("tables.json"))
r = tskit_sample_nodes_by_ploidy(tables, ploidy=2)
print(r["groups"])
```

- **Example output:**

```text
[[0, 1], [2, 3]]
```

## 10) `dnakit.population.tskit_ibd_segments()`

- **Purpose:** Identity-by-descent segments between sample pairs from the genealogy.
- **Method:** Calls tskit `TreeSequence.ibd_segments`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_ibd_segments(tables, *, within=None, between=None, max_time=None, min_span=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `within` or `between`, `max_time`, `min_span`.
- **Output:** `segments` rows with `node_a`, `node_b`, `left`, `right`, `mrca`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_ibd_segments

tables = json.load(open("tables.json"))
r = tskit_ibd_segments(tables, within=[0, 1, 2])
print(r["num_segments"])
for s in r["segments"]:
    print(s)
```

- **Example output:**

```text
6
{'node_a': 0, 'node_b': 1, 'left': 0.0, 'right': 50.0, 'mrca': 4}
{'node_a': 0, 'node_b': 1, 'left': 50.0, 'right': 100.0, 'mrca': 9}
{'node_a': 0, 'node_b': 2, 'left': 50.0, 'right': 100.0, 'mrca': 7}
{'node_a': 0, 'node_b': 2, 'left': 0.0, 'right': 50.0, 'mrca': 6}
{'node_a': 1, 'node_b': 2, 'left': 0.0, 'right': 50.0, 'mrca': 6}
{'node_a': 1, 'node_b': 2, 'left': 50.0, 'right': 100.0, 'mrca': 9}
```

## 11) `dnakit.population.tskit_read_hapmap()`

- **Purpose:** Read a HapMap-style recombination map file.
- **Method:** Calls tskit `RateMap.read_hapmap`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_read_hapmap(source, *, sequence_length=None, has_header=True, position_col=None, rate_col=None, map_col=None)`.
- **Input:** `source` file, `sequence_length`, `has_header`, `position_col`, `rate_col`, `map_col`.
- **Output:** `intervals` with `left`, `right`, `rate`, plus `total_mass` and `mean_rate`.
- **Example code:**

```python
from dnakit.population.tskit_p31_functions import tskit_read_hapmap

r = tskit_read_hapmap("map.txt")
print(r["sequence_length"], r["total_mass"], r["mean_rate"])
for i in r["intervals"]:
    print(i)
```

- **Example output:**

```text
2000000.0 0.03 1.5e-08
{'left': 0.0, 'right': 1000000.0, 'rate': 1e-08}
{'left': 1000000.0, 'right': 2000000.0, 'rate': 1.9999999999999997e-08}
```

## 12) `dnakit.population.tskit_recombination_map_lookup()`

- **Purpose:** Look up recombination rate or cumulative genetic position at positions.
- **Method:** Calls tskit `RateMap.get_rate` / `get_cumulative_mass`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_recombination_map_lookup(positions, *, quantity='rate', hapmap_file=None, map_positions=None, map_rates=None, sequence_length=None, has_header=True)`.
- **Input:** `positions`, `quantity`, and the map as `hapmap_file` or `map_positions` with `map_rates`.
- **Output:** `values`, one per position.
- **Example code:**

```python
from dnakit.population.tskit_p31_functions import tskit_recombination_map_lookup

r = tskit_recombination_map_lookup([500000, 1500000], hapmap_file="map.txt")
print(r["values"])
r = tskit_recombination_map_lookup([1000000], quantity="cumulative_mass", hapmap_file="map.txt")
print(r["values"])
```

- **Example output:**

```text
[1e-08, 1.9999999999999997e-08]
[0.01]
```
