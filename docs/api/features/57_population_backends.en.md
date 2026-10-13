# Optional tskit and dadi population-genetics backends

This group of functions in `dnakit.population` is a set of adapters over optional packages. [tskit](https://tskit.dev/tskit) builds tree sequences from explicit tables, reads and writes `.trees` files, and computes per-site diversity, Tajima's D, Fst, and the allele frequency spectrum. [dadi](https://github.com/RyanGutenkunst/dadi) provides the expected spectrum for a constant-size neutral model, spectrum folding, and the Poisson log-likelihood of an observed spectrum. Sequence simulation and fitting demographic history models are not included.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command.

## Example data

The tskit examples in items 1 to 7 read the table below. It has 4 samples, 1 tree, and 3 sites. It only demonstrates the interface, and its statistics do not describe a real population. The dadi examples in items 8 to 11 use plain lists and need no file.

**Example file:** `tables.json`

```json
{
  "sequence_length": 100,
  "nodes": [{"time": 0, "sample": true}, {"time": 0, "sample": true}, {"time": 0, "sample": true}, {"time": 0, "sample": true}, {"time": 1}, {"time": 2}, {"time": 3}],
  "edges": [{"left": 0, "right": 100, "parent": 4, "child": 0}, {"left": 0, "right": 100, "parent": 4, "child": 1}, {"left": 0, "right": 100, "parent": 5, "child": 2}, {"left": 0, "right": 100, "parent": 5, "child": 3}, {"left": 0, "right": 100, "parent": 6, "child": 4}, {"left": 0, "right": 100, "parent": 6, "child": 5}],
  "sites": [{"position": 10, "ancestral_state": "A"}, {"position": 20, "ancestral_state": "A"}, {"position": 30, "ancestral_state": "A"}],
  "mutations": [{"site": 0, "node": 0, "derived_state": "G"}, {"site": 1, "node": 4, "derived_state": "G"}, {"site": 2, "node": 5, "derived_state": "G"}]
}
```
## 1) `dnakit.population.build_tree_sequence()`

- **Purpose:** Build a tskit tree sequence from explicit node, edge, site, and mutation rows, and return a summary.
- **Method:** Write the rows of each table into a tskit `TableCollection`, call `sort()`, then call `tree_sequence()`. Intervals are half-open `[left, right)` in the same units as `sequence_length`. Invalid tables raise an error from tskit.
- **API:** `dnakit.population.build_tree_sequence(tables, *, max_rows=100000)`.
- **Input:** `tables` is a mapping with `sequence_length` (greater than 0), `nodes` (each row has `time` and optionally `sample`), `edges` (`left`, `right`, `parent`, `child`), `sites` (`position`, `ancestral_state`), and `mutations` (`site`, `node`, `derived_state`, and optionally `parent` and `time`). `max_rows` limits the number of rows in each table.
- **Output:** A dictionary with `num_nodes`, `num_edges`, `num_samples`, `num_trees`, `num_sites`, `num_mutations`, `num_populations`, `sequence_length`, and `provenance`.
- **Example code:**

```python
import json
from dnakit.population import build_tree_sequence

tables = json.load(open("tables.json"))
result = build_tree_sequence(tables)
print(result["num_trees"], result["num_samples"], result["num_sites"])
```

- **Example output:**

```text
1 4 3
```

## 2) `dnakit.population.write_tree_sequence()`

- **Purpose:** Build a tree sequence from tables and write it as a tskit `.trees` file.
- **Method:** Use the same building and validation rules as `build_tree_sequence()`, then call `TreeSequence.dump()`. An existing file is never overwritten, and the parent directory must already exist.
- **API:** `dnakit.population.write_tree_sequence(tables, destination, *, max_rows=100000)`.
- **Input:** The same `tables` as `build_tree_sequence()`; `destination` is the output path.
- **Output:** The written path, the file size in bytes, and the same summary fields as `build_tree_sequence()`.
- **Example code:**

```python
import json
from pathlib import Path
from dnakit.population import write_tree_sequence

Path("sample.trees").unlink(missing_ok=True)  # the function never overwrites an existing file
tables = json.load(open("tables.json"))
result = write_tree_sequence(tables, "sample.trees")
print(result["destination"].endswith("sample.trees"), result["num_mutations"])
```

- **Example output:**

```text
True 3
```

## 3) `dnakit.population.read_tree_sequence()`

- **Purpose:** Read a tskit `.trees` file and return a summary.
- **Method:** Check that the file size does not exceed `max_bytes`, then call `tskit.load()`. A file that is not in tskit format raises an error.
- **API:** `dnakit.population.read_tree_sequence(source, *, max_bytes=100000000)`.
- **Input:** The path to a `.trees` file; `max_bytes` is the size limit.
- **Output:** A dictionary with `source`, `bytes`, and the same summary fields as `build_tree_sequence()`.
- **Example code:**

```python
from dnakit.population import read_tree_sequence

result = read_tree_sequence("sample.trees")
print(result["num_trees"], result["num_edges"])
```

- **Example output:**

```text
1 6
```

## 4) `dnakit.population.calculate_diversity()`

- **Purpose:** Compute nucleotide diversity per base (pi) for each sample set.
- **Method:** Call tskit's `diversity()` in site mode, and divide the result by `sequence_length`. When `sample_sets` is omitted, all samples form one set named `"all"`.
- **API:** `dnakit.population.calculate_diversity(tables, *, sample_sets=None, max_rows=100000)`.
- **Input:** `tables`; `sample_sets` is `{name: [sample node IDs, ...]}`. Each set needs at least two distinct sample nodes.
- **Output:** `diversity_per_base` (one value per set), `sequence_length`, and `provenance`.
- **Example code:**

```python
import json
from dnakit.population import calculate_diversity

tables = json.load(open("tables.json"))
result = calculate_diversity(tables, sample_sets={"A": [0, 1], "B": [2, 3]})
print({name: round(value, 6) for name, value in result["diversity_per_base"].items()})
```

- **Example output:**

```text
{'A': 0.01, 'B': 0.0}
```

## 5) `dnakit.population.calculate_tajimas_d()`

- **Purpose:** Compute Tajima's D and the number of segregating sites for each sample set.
- **Method:** Call tskit's `Tajimas_D()` and `segregating_sites()`. When tskit returns NaN (for example, when there are no segregating sites), the value is reported as `None`.
- **API:** `dnakit.population.calculate_tajimas_d(tables, *, sample_sets=None, max_rows=100000)`.
- **Input:** The same as `calculate_diversity()`.
- **Output:** `tajimas_d`, `segregating_sites`, and `provenance`.
- **Example code:**

```python
import json
from dnakit.population import calculate_tajimas_d

tables = json.load(open("tables.json"))
result = calculate_tajimas_d(tables)
print(result["segregating_sites"]["all"], round(result["tajimas_d"]["all"], 4))
```

- **Example output:**

```text
3 1.0898
```

## 6) `dnakit.population.calculate_fst()`

- **Purpose:** Compute Fst between two sample sets.
- **Method:** Call tskit's `Fst()` in site mode. tskit's estimator is `1 - 2(d_X + d_Y) / (d_X + 2 d_XY + d_Y)`, where `d_X` and `d_Y` are within-set diversities and `d_XY` is the between-set divergence. When the denominator is zero, the value is `None`.
- **API:** `dnakit.population.calculate_fst(tables, *, sample_sets, max_rows=100000)`.
- **Input:** `tables`; `sample_sets` must contain exactly two sets.
- **Output:** `fst`, the name and sample count of each population, and `provenance`.
- **Example code:**

```python
import json
from dnakit.population import calculate_fst

tables = json.load(open("tables.json"))
result = calculate_fst(tables, sample_sets={"A": [0, 1], "B": [2, 3]})
print(round(result["fst"], 4))
```

- **Example output:**

```text
0.6667
```

## 7) `dnakit.population.calculate_site_frequency_spectrum()`

- **Purpose:** Count the site frequency spectrum for each sample set.
- **Method:** Call tskit's `allele_frequency_spectrum()`, returning counts rather than frequencies. `spectrum[i]` is the number of sites where the derived allele is carried by `i` samples in the set. With `folded=True`, the minor-allele spectrum is returned, so the counts for `i` and `n - i` are merged into `min(i, n - i)`.
- **API:** `dnakit.population.calculate_site_frequency_spectrum(tables, *, sample_sets=None, folded=False, max_rows=100000)`.
- **Input:** `tables`; `sample_sets` as in `calculate_diversity()`; `folded` is `True` or `False`.
- **Output:** `spectra` (one list of length `n + 1` per set), `folded`, and `provenance`.
- **Example code:**

```python
import json
from dnakit.population import calculate_site_frequency_spectrum

tables = json.load(open("tables.json"))
result = calculate_site_frequency_spectrum(tables)
print(result["spectra"]["all"])
```

- **Example output:**

```text
[0, 1, 2, 0, 0]
```

## 8) `dnakit.population.expected_neutral_sfs()`

- **Purpose:** Compute the expected site frequency spectrum under the standard neutral model for a constant-size population.
- **Method:** Call dadi's `Demographics1D.snm()`, which integrates numerically on a grid of `grid_points` points. The expected values approximate `theta / i`, and the result is multiplied by `theta`. With `folded=True`, dadi's `fold()` is applied. Classes 0 and n are monomorphic, and dadi masks them, so they are returned as `None`. A denser grid reduces the error.
- **API:** `dnakit.population.expected_neutral_sfs(sample_size, *, theta=1.0, folded=False, grid_points=1000, max_sample_size=1000)`.
- **Input:** The number of sampled chromosomes `sample_size` (at least 2); `theta` is a number greater than zero; `grid_points` is an integer from 50 to 10000.
- **Output:** `spectrum` (a list of length `sample_size + 1`, with masked entries as `None`), `theta`, `folded`, and `provenance`.
- **Example code:**

```python
from dnakit.population import expected_neutral_sfs

result = expected_neutral_sfs(4, theta=1.0)
print([None if value is None else round(value, 3) for value in result["spectrum"]])
```

- **Example output:**

```text
[None, 1.0, 0.5, 0.333, None]
```

## 9) `dnakit.population.fold_sfs()`

- **Purpose:** Fold an unfolded site frequency spectrum into the minor-allele spectrum.
- **Method:** Use dadi's `Spectrum(...).fold()`. Classes 0 and n are masked by dadi and do not take part in the fold. Entries that are masked after folding are returned as `None`.
- **API:** `dnakit.population.fold_sfs(spectrum, *, max_sample_size=1000)`.
- **Input:** A list of non-negative counts of length `n + 1`, where `n` is at least 2.
- **Output:** `sample_size`, `spectrum` (the folded list, with masked entries as `None`), and `provenance`.
- **Example code:**

```python
from dnakit.population import fold_sfs

result = fold_sfs([0, 2, 1, 1, 0])
print(result["spectrum"])
```

- **Example output:**

```text
[None, 3.0, 1.0, None, None]
```

## 10) `dnakit.population.estimate_sfs_theta()`

- **Purpose:** Estimate the maximum-likelihood `theta` of an observed spectrum under the neutral model.
- **Method:** Use dadi's `Inference.optimal_sfs_scaling()` on classes 1 to n-1. In theory, this equals the observed total divided by the model total.
- **API:** `dnakit.population.estimate_sfs_theta(observed, *, grid_points=1000, max_sample_size=1000)`.
- **Input:** Observed counts `observed` of length `n + 1`. The sum over classes 1 to n-1 must be greater than zero.
- **Output:** `sample_size`, `theta`, and `provenance`.
- **Example code:**

```python
from dnakit.population import estimate_sfs_theta

result = estimate_sfs_theta([0, 2, 1, 1, 0])
print(round(result["theta"], 3))
```

- **Example output:**

```text
2.182
```

## 11) `dnakit.population.calculate_sfs_log_likelihood()`

- **Purpose:** Compute the Poisson log-likelihood of an observed spectrum under the neutral model.
- **Method:** Use dadi's `Inference.ll()` on classes 1 to n-1, including the `log(x!)` term. When `theta` is omitted, the maximum-likelihood `theta` is used first, which gives the same result as dadi's `ll_multinom()`.
- **API:** `dnakit.population.calculate_sfs_log_likelihood(observed, *, theta=None, grid_points=1000, max_sample_size=1000)`.
- **Input:** `observed` as in `estimate_sfs_theta()`; `theta` may be omitted, or given as a number greater than zero.
- **Output:** `theta`, `theta_source` (`"maximum likelihood"` or `"given"`), `log_likelihood`, and `provenance`.
- **Example code:**

```python
from dnakit.population import calculate_sfs_log_likelihood

result = calculate_sfs_log_likelihood([0, 2, 1, 1, 0])
print(result["theta_source"], round(result["log_likelihood"], 3))
```

- **Example output:**

```text
maximum likelihood -3.364
```

## Installing the optional dependencies

```bash
python -m pip install "dnakit[bio-tskit]"    # tskit: build/read/write_tree_sequence and calculate_* population statistics
python -m pip install "dnakit[bio-dadi]"     # dadi: expected_neutral_sfs and the other spectrum functions
```

Project pages: [tskit](https://tskit.dev/tskit), [dadi](https://github.com/RyanGutenkunst/dadi).

## How these functions are verified

The test file is `others/tests/unit/test_population_backends.py`, with 42 tests. They check the following:

- The tree-sequence summary, diversity, Fst, and the allele frequency spectrum match direct calls to tskit. Diversity, Tajima's D, and Fst also match values calculated by hand from the formulas. Tajima's D is additionally compared with an exact fraction-based calculation of Tajima (1989).
- Reversing the order of the site rows does not change the statistics after tskit sorts the tables. A file written by `write_tree_sequence()` can be read back by `read_tree_sequence()`.
- The dadi neutral expected spectrum matches `theta / i` within a relative error of 1e-4, and it matches a direct call to `Demographics1D.snm()`. The folded spectrum matches a hand fold. The maximum-likelihood `theta` for the example matches `24/11`. The log-likelihood matches the Poisson formula and `Inference.ll_multinom()`.
- Booleans, out-of-range values, unknown keys, non-sample nodes, values above `max_rows` or `max_sample_size`, and existing output files are all rejected.
- A missing `tskit` or `dadi` package raises an error that includes the matching extra name.

## Current scope

This batch covers tskit tree-sequence reading and writing, per-site diversity, Tajima's D, Fst between two populations, and the site frequency spectrum, together with dadi's constant-size neutral model. Tree sequences come only from explicit tables or an existing `.trees` file. Simulation is not included (msprime is not wrapped).

Not wrapped from tskit: branch-mode and windowed statistics, divergence, linkage disequilibrium, genotype matrices, tree traversal, and VCF export. Not wrapped from dadi: demographic models, reading frequency-spectrum data dictionaries, parameter optimization, Godambe variance, and two-dimensional spectra. The dadi expected spectrum is a numerical integration, and its error decreases as `grid_points` increases.
