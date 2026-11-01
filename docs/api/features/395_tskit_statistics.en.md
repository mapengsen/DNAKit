# tskit sample-set, relatedness, trait and coalescence statistics

These functions in `dnakit.population` wrap the statistics of [tskit](https://github.com/tskit-dev/tskit): linkage disequilibrium, Y1/Y2/Y3/f4, segregating sites, general statistics, genetic relatedness, trait association, mean descendants, genealogical nearest neighbours and pair coalescence.

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

## 1) `dnakit.population.tskit_ld_r2()`

- **Purpose:** Linkage disequilibrium r2 between sites as a matrix, a pair or an array.
- **Method:** Calls tskit `tskit.LdCalculator.r2_matrix` / `r2` / `r2_array`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_ld_r2(tables, *, output='matrix', site_a=None, site_b=None, focal_site=None, direction=1, max_mutations=None, max_distance=None, max_sites=None)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `output`, `site_a`, `site_b`, `focal_site`, `direction`, `max_mutations`, `max_distance`, `max_sites`.
- **Output:** `values` (matrix, number or array).
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_ld_r2

tables = json.load(open("tables.json"))
print(tskit_ld_r2(tables, output="pair", site_a=0, site_b=1)["values"])
print(tskit_ld_r2(tables, output="array", focal_site=0)["values"])
print(tskit_ld_r2(tables)["values"][0])
```

- **Example output:**

```text
1.0
[1.0, 0.0, 0.0]
[1.0, 1.0, 0.0, 0.0]
```

## 2) `dnakit.population.tskit_sample_set_statistic()`

- **Purpose:** The multi-set statistics Y1, Y2, Y3 and f4.
- **Method:** Calls tskit `TreeSequence.Y1` / `Y2` / `Y3` / `f4`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_sample_set_statistic(tables, *, statistic, sample_sets=None, indexes=None, windows=None, mode='site', span_normalise=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `statistic`, `sample_sets`, `indexes`, `windows`, `mode`, `span_normalise`.
- **Output:** `values` (windows by sets or index tuples) plus `mode`, `windows`, `span_normalise`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_sample_set_statistic

tables = json.load(open("tables.json"))
sets = [[0], [1], [2], [3]]
print(tskit_sample_set_statistic(tables, statistic="f4", sample_sets=sets, indexes=[[0, 1, 2, 3]])["values"])
print(tskit_sample_set_statistic(tables, statistic="Y3", sample_sets=sets, indexes=[[0, 1, 2]])["values"])
print(tskit_sample_set_statistic(tables, statistic="Y1", sample_sets=[[0, 1, 2]])["values"])
```

- **Example output:**

```text
[0.02]
[0.0]
[0.013333333333333332]
```

## 3) `dnakit.population.tskit_segregating_sites()`

- **Purpose:** Segregating-site density (or count) per window and sample set.
- **Method:** Calls tskit `TreeSequence.segregating_sites`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_segregating_sites(tables, *, sample_sets=None, windows=None, mode='site', span_normalise=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `sample_sets`, `windows`, `mode`, `span_normalise`.
- **Output:** `values` (windows by sets or index tuples) plus `mode`, `windows`, `span_normalise`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_segregating_sites

tables = json.load(open("tables.json"))
print(tskit_segregating_sites(tables, sample_sets=[[0, 1], [2, 3]], span_normalise=False)["values"])
print(tskit_segregating_sites(tables, windows=[0, 50, 100], span_normalise=False)["values"])
```

- **Example output:**

```text
[2.0, 2.0]
[[2.0], [2.0]]
```

## 4) `dnakit.population.tskit_general_stat()`

- **Purpose:** A windowed statistic from sample weights or sample-set counts, with a named summary function.
- **Method:** Calls tskit `TreeSequence.general_stat` / `sample_count_stat`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_general_stat(tables, *, weights=None, sample_sets=None, summary='identity', polarised=False, strict=True, windows=None, mode='site', span_normalise=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `weights` or `sample_sets`, `summary` (identity, square, product), `polarised`, `strict`, `windows`, `mode`, `span_normalise`.
- **Output:** `values` (windows by sets or index tuples) plus `mode`, `windows`, `span_normalise`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_general_stat

tables = json.load(open("tables.json"))
w = [[1.0], [-1.0], [2.0], [-2.0]]
print(tskit_general_stat(tables, weights=w, summary="square", polarised=True)["values"])
print(tskit_general_stat(tables, sample_sets=[[0, 1], [2, 3]], summary="product", polarised=True, strict=False)["values"])
```

- **Example output:**

```text
[0.18]
[0.02]
```

## 5) `dnakit.population.tskit_genetic_relatedness()`

- **Purpose:** Genetic relatedness between sample sets, as pairs or a full matrix.
- **Method:** Calls tskit `TreeSequence.genetic_relatedness` / `genetic_relatedness_matrix`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_genetic_relatedness(tables, *, sample_sets=None, indexes=None, windows=None, mode='site', span_normalise=True, centre=True, proportion=True, polarised=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `sample_sets`, `indexes`, `centre`, `proportion`, `polarised`, `windows`, `mode`, `span_normalise`.
- **Output:** `values` (windows by sets or index tuples) plus `mode`, `windows`, `span_normalise`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genetic_relatedness

tables = json.load(open("tables.json"))
sets = [[0], [1], [2], [3]]
print(tskit_genetic_relatedness(tables, sample_sets=sets, indexes=[[0, 1], [0, 2]])["values"])
print(tskit_genetic_relatedness(tables, sample_sets=sets)["values"])
```

- **Example output:**

```text
[0.0, 0.0]
[[0.01, -0.0, -0.0, -0.009999999999999998], [-0.0, 0.01, -0.009999999999999998, -0.0], [-0.0, -0.009999999999999998, 0.01, -0.0], [-0.009999999999999998, -0.0, -0.0, 0.01]]
```

## 6) `dnakit.population.tskit_genetic_relatedness_vector()`

- **Purpose:** Product of the genetic relatedness matrix with a weight matrix (branch mode).
- **Method:** Calls tskit `TreeSequence.genetic_relatedness_vector`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_genetic_relatedness_vector(tables, *, weights, windows=None, span_normalise=True, centre=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `weights` (one row per sample), `centre`, `windows`, `span_normalise`.
- **Output:** `values`, same shape as `weights` per window.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genetic_relatedness_vector

tables = json.load(open("tables.json"))
r = tskit_genetic_relatedness_vector(tables, weights=[[1.0], [2.0], [0.5], [3.0]])
print(r["values"])
```

- **Example output:**

```text
[[-4.03125], [2.59375], [-4.46875], [5.90625]]
```

## 7) `dnakit.population.tskit_genetic_relatedness_weighted()`

- **Purpose:** Genetic relatedness between pairs of weight columns.
- **Method:** Calls tskit `TreeSequence.genetic_relatedness_weighted`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_genetic_relatedness_weighted(tables, *, weights, indexes=None, windows=None, mode='site', span_normalise=True, polarised=False, centre=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `weights`, `indexes`, `polarised`, `centre`, `windows`, `mode`, `span_normalise`.
- **Output:** `values` (windows by sets or index tuples) plus `mode`, `windows`, `span_normalise`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genetic_relatedness_weighted

tables = json.load(open("tables.json"))
r = tskit_genetic_relatedness_weighted(tables, weights=[[1.0, 0.5], [2.0, 1.0], [0.5, 3.0], [3.0, 2.0]], indexes=[[0, 1]])
print(r["values"])
```

- **Example output:**

```text
[0.0]
```

## 8) `dnakit.population.tskit_trait_statistic()`

- **Purpose:** Association of genotype with traits: correlation, covariance or linear model.
- **Method:** Calls tskit `TreeSequence.trait_correlation` / `trait_covariance` / `trait_linear_model`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_trait_statistic(tables, *, statistic, weights, covariates=None, windows=None, mode='site', span_normalise=True)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `statistic`, `weights` (sample by trait), `covariates`, `windows`, `mode`, `span_normalise`.
- **Output:** `values` (windows by sets or index tuples) plus `mode`, `windows`, `span_normalise`.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_trait_statistic

tables = json.load(open("tables.json"))
w = [[1.0], [2.0], [0.5], [3.0]]
for stat in ("correlation", "covariance", "linear_model"):
    print(stat, tskit_trait_statistic(tables, statistic=stat, weights=w)["values"])
```

- **Example output:**

```text
correlation [0.01694915254237288]
covariance [0.006944444444444444]
linear_model [0.0625]
```

## 9) `dnakit.population.tskit_mean_descendants()`

- **Purpose:** Mean number of descendant samples of each sample set for every node.
- **Method:** Calls tskit `TreeSequence.mean_descendants`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_mean_descendants(tables, *, sample_sets)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `sample_sets`.
- **Output:** `values`, one row per node and one column per sample set.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_mean_descendants

tables = json.load(open("tables.json"))
r = tskit_mean_descendants(tables, sample_sets=[[0, 1], [2, 3]])
for row in r["values"]:
    print(row)
```

- **Example output:**

```text
[1.0, 0.0]
[1.0, 0.0]
[0.0, 1.0]
[0.0, 1.0]
[2.0, 0.0]
[0.0, 2.0]
[2.0, 2.0]
[1.0, 1.0]
[1.0, 1.0]
[2.0, 2.0]
```

## 10) `dnakit.population.tskit_genealogical_nearest_neighbours()`

- **Purpose:** Share of each focal node's genealogical nearest neighbours that lie in each sample set.
- **Method:** Calls tskit `TreeSequence.genealogical_nearest_neighbours`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_genealogical_nearest_neighbours(tables, *, focal, sample_sets)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `focal`, `sample_sets` (non-overlapping).
- **Output:** `values`, one row per focal node.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_genealogical_nearest_neighbours

tables = json.load(open("tables.json"))
r = tskit_genealogical_nearest_neighbours(tables, focal=[0, 3], sample_sets=[[1], [2]])
print(r["values"])
```

- **Example output:**

```text
[[0.5, 0.5], [0.5, 0.5]]
```

## 11) `dnakit.population.tskit_pair_coalescence()`

- **Purpose:** Pair coalescence counts, quantiles or rates over time.
- **Method:** Calls tskit `TreeSequence.pair_coalescence_counts` / `pair_coalescence_quantiles` / `pair_coalescence_rates`; the result is converted to plain Python values.
- **API:** `dnakit.population.tskit_pair_coalescence(tables, *, output='rates', time_windows='nodes', quantiles=None, sample_sets=None, indexes=None, windows=None, span_normalise=True, pair_normalise=False)`.
- **Input:** `tables` (explicit rows or a `.trees` path), `output`, `time_windows`, `quantiles`, `sample_sets`, `indexes`, `windows`, `span_normalise`, `pair_normalise`.
- **Output:** `values` and the `time_windows` or `quantiles` used.
- **Example code:**

```python
import json
from dnakit.population.tskit_p31_functions import tskit_pair_coalescence

tables = json.load(open("tables.json"))
sets = [[0, 1, 2, 3]]
print(tskit_pair_coalescence(tables, output="counts", sample_sets=sets)["values"])
print(tskit_pair_coalescence(tables, output="rates", time_windows=[0, 2], sample_sets=sets)["values"])
print(tskit_pair_coalescence(tables, output="quantiles", quantiles=[0.5, 1.0], sample_sets=sets)["values"])
```

- **Example output:**

```text
[0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 2.0, 0.5, 0.5, 2.0]
[0.09116077839697735, 0.4761904761904763]
[4.0, 5.0]
```
