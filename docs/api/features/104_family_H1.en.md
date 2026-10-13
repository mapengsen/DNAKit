# H1 family-variant review: dadi, tskit and scikit-allel

This page records the review of the `family-variant` ledger rows in group H1. Each of these rows had been classed as an alias, a re-export or a parameter-selected branch of a function that is already wrapped. Every claim was checked, and the functions that have their own DNA-level meaning and output were wrapped. The new functions live in `src/dnakit/population/dadi_family_functions.py`, `src/dnakit/population/tskit_family_functions.py` and `src/dnakit/variants/allel_family_functions.py`. All of them are exported from `dnakit.population` and `dnakit.variants`, and the examples import from those packages.

- 26 functions are new: 6 for dadi, 2 for tskit and 18 for scikit-allel. Each has an example with real output below.
- The tests are in `others/tests/unit/test_family_H1.py`. Each new function is checked against the upstream call or a hand calculation. For a sample of the rows kept as `family-variant`, the existing wrapper is called with the claimed parameter and compared with the upstream variant.
- The example outputs come from running the code, so the documentation checker should report 0 mismatches.

## Review results

| Package | family-variant before | duplicate-of | family-variant (verified) | wrapped | Other statuses |
| --- | --- | --- | --- | --- | --- |
| dadi | 121 | 1 | 108 | 10 | broken-upstream 2 |
| tskit | 11 | 0 | 9 | 2 | - |
| allel | 28 | 1 | 0 | 25 | broken-upstream 1; needs-binary 1 |

In total 160 rows were checked: 117 stay `family-variant`, 37 became `wrapped`, 2 became `duplicate-of`, and the other 4 became `broken-upstream` or `needs-binary`.

### Rows that stay family-variant (verified)

A row stays only when the claimed parameter is in the existing wrapper's signature (checked with `inspect.signature`) and the wrapper calls exactly the upstream object of the row.

| Wrapper and parameter | Rows | Upstream functions |
| --- | --- | --- |
| `expected_sfs_1d(model=...)` | 5 | `Demographics1D.growth`, `bottlegrowth`, `bottlegrowth_1d`, `three_epoch`, `three_epoch_inbreeding` |
| `expected_sfs_2d(model=...)` | 42 | the other `Demographics2D` models (`no_mig`, `IM`, `founder_asym`, ...) |
| `expected_sfs_3d(model=...)` | 33 | the other `Demographics3D` models (`out_of_africa`, ...) |
| `expected_sfs_selection(model=...)` | 20 | the other `DFE.DemogSelModels` models |
| `ms_core_command_2d(model='IM' / 'IM_pre')` | 2 | `Demographics2D.IM_mscore`, `IM_pre_mscore` |
| `nested_model_statistic(test='wald' / 'score')` | 2 | `Godambe.Wald_stat`, `score_stat` |
| `fit_demographic_model(method='nelder-mead' / 'powell')` | 2 | `Inference.optimize_log_fmin`, `optimize_log_powell` |
| `sfs_residuals(residual='poisson')` | 1 | `Inference.linear_Poisson_residual` |
| `read_vcf_data(parser='cyvcf2')` | 1 | `Misc.make_data_dict_vcf_cyvcf2` |
| `read_text_table(table=...)` | 6 | `tskit.parse_individuals`, `parse_edges`, `parse_sites`, `parse_mutations`, `parse_populations`, `parse_migrations` |
| `random_tree_newick(shape='star' / 'comb' / 'balanced')` | 3 | `tskit.combinatorics.generate_star`, `generate_comb`, `generate_balanced` |

For the 1D, 2D, 3D and selection models of dadi (100 rows) every model name was run through the wrapper: it is in the wrapper's fixed list, `getattr` returns the object of the row, and the call succeeds on a small sample size and grid.

### Claims that were wrong

- Five dadi likelihood functions (`ll_multinom`, `minus_ll`, `minus_ll_multinom`, `ll_per_bin`, `ll_multinom_per_bin`) and `optimally_scaled_sfs` were called variants of `calculate_sfs_log_likelihood` and `estimate_sfs_theta`, and the ledger itself said “not exposed”. Those two wrappers work on the neutral model and a one-dimensional spectrum only, and no parameter reaches the functions, so they became `wrapped` (`calculate_model_log_likelihood`, `scale_model_sfs_to_data`).
- dadi `optimize`, `optimize_cons`, `optimize_grid` and `opt` were called variants of `fit_demographic_model`. Its `method` has four log-scale values only, and no constraint, grid or NLopt option, so these four rows became `wrapped`. `NLopt_mod.opt` is the same object as `Inference.opt` (`is` is true) and became `duplicate-of`.
- `optimize_log_lbfgsb` was said to be reached by `fit_demographic_model(method='lbfgsb')`, but dadi 2.4.4 passes `iprint` to `scipy.optimize.fmin_l_bfgs_b`, which SciPy 1.18.1 no longer accepts, so upstream raises TypeError (`optimize_lbfgsb` fails the same way). Both rows became `broken-upstream`. Note that the docstring of `fit_demographic_model` still lists `lbfgsb`, and in this environment it always ends in a ConfigurationError.
- Of the 28 scikit-allel rows, only the `average_patterson_d` claim holds, and it is not a parameter variant: it is the same object as `blockwise_patterson_d` (`duplicate-of`). For all the others, the windowed and per-variant forms, joint spectra, folded and scaled forms, jackknife forms, randomized PCA, phasing steps and field-table forms, no parameter of an existing wrapper reaches them, and 25 rows became `wrapped`. For example `calculate_fst` has the single parameter `method` and returns only an average.
- `allel.roh_poissonhmm` needs the optional package protopunica, which is not in the full environment and cannot be installed here, so it became `needs-binary`. `allel.write_vcf` is documented upstream as preliminary, and a check showed that it writes malformed VCF (ALT as `G,`, a missing QUAL as `nan`, every FILTER as PASS), so it became `broken-upstream`.
- tskit `write_nexus` was said to be reached by `export_sample_sequences(file_format='nexus')`, but that wrapper fixes `include_trees=False` and cannot write the TREES block. `tree_count_topologies` was called the single-tree version of `count_sample_topologies`, which has no `tree_index` and evaluates every tree. Both rows became `wrapped`.

### Other findings (outside the rows of this group)

- `allel.fold_sfs` gives a wrong result for an even number of chromosomes: for `[1, 2, 7, 14, 11, 5, 0]` (n=6) it returns `[6, 13, 21]`, while the correct folded spectrum is `[1, 7, 18, 14]` (class i joined with class n-i, which `sfs_folded` gives, and dadi's fold gives 7, 18, 14 for the unmasked classes 1 to 3). The ledger row that marks it `duplicate-of` the dadi fold does not hold for even n.
- Existing wrappers such as `expected_sfs_selection` convert only ValueError, TypeError, IndexError and KeyError to ConfigurationError. dadi's `IM_pre_sel` raises ZeroDivisionError for some parameters, and it reaches the caller unchanged.

## Shared conventions

- The `verbose` argument of every public function is keyword-only and defaults to `True`.
- Return values are JSON-serializable dicts with `status`, `operation`, the result fields and `provenance`.
- A statistic that cannot be computed is `null`, and so are the entries dadi masks.
- dadi, tskit and scikit-allel are imported only when a function runs. A missing dependency raises `BackendUnavailableError`.
- Invalid input raises `ConfigurationError`. Output files never overwrite an existing file.
- Randomized PCA fixes its randomness through `seed` only and changes no global random state.

## 1) `calculate_model_log_likelihood`

- **Purpose:** Log-likelihood of observed counts given any model spectrum.
- **Method:** Builds dadi spectra from `model` and `data` (corners masked) and calls `Inference.ll`. `multinom=True` uses `ll_multinom` (the model is first rescaled to the optimal theta), `negate=True` uses `minus_ll` or `minus_ll_multinom`, and `per_bin=True` uses `ll_per_bin` or `ll_multinom_per_bin` and returns the value of every entry. `negate` and `per_bin` cannot be combined. The existing `calculate_sfs_log_likelihood` only gives the Poisson total for the neutral model and has none of these options.
- **API:** `calculate_model_log_likelihood(model[required], data[required], *, multinom[optional], per_bin[optional], negate[optional], verbose[optional])`
- **Input:** `model` and `data` are nested lists of the same shape (one to three populations, each axis at least 3 long). `data` must be non-negative.
- **Output:** A dict with `dadi_function` (the dadi function that ran), `multinom`, `per_bin`, `negate`, `shape`, and either `log_likelihood` (a number) or `log_likelihood_per_bin` (a nested list with null for masked entries).

**Example code:**

```python
from dnakit.population import calculate_model_log_likelihood

model = [0.0, 10.0, 5.0, 3.3, 2.5, 2.0, 1.7, 1.4, 1.2, 1.0, 0.0]
data = [0.0, 12.0, 4.0, 4.0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 0.0]
plain = calculate_model_log_likelihood(model, data, verbose=False)
multi = calculate_model_log_likelihood(model, data, multinom=True, verbose=False)
negated = calculate_model_log_likelihood(model, data, negate=True, verbose=False)
per_bin = calculate_model_log_likelihood(model, data, per_bin=True, verbose=False)
print(plain["dadi_function"], round(plain["log_likelihood"], 4))
print(multi["dadi_function"], round(multi["log_likelihood"], 4))
print(negated["dadi_function"], round(negated["log_likelihood"], 4))
print([None if v is None else round(v, 3) for v in per_bin["log_likelihood_per_bin"]])
```

**Example output:**

```text
ll -13.0735
ll_multinom -13.0593
minus_ll 13.0735
[None, -2.356, -1.74, -1.702, -1.361, -1.307, -1.169, -1.42, -1.018, -1.0, None]
```

## 2) `scale_model_sfs_to_data`

- **Purpose:** Scale a model spectrum by the maximum-likelihood theta for the data.
- **Method:** Calls `Inference.optimally_scaled_sfs` for the scaled spectrum and `Inference.optimal_sfs_scaling` for the factor theta (the data total over the model total, on unmasked entries). The existing `estimate_sfs_theta` returns theta for the neutral model only and not the scaled spectrum.
- **API:** `scale_model_sfs_to_data(model[required], data[required], *, verbose[optional])`
- **Input:** `model` and `data` as above, with the same shape.
- **Output:** A dict with `shape`, `theta` and `scaled_spectrum` (a nested list with null for masked entries).

**Example code:**

```python
from dnakit.population import scale_model_sfs_to_data

model = [0.0, 10.0, 5.0, 3.3, 2.5, 2.0, 1.7, 1.4, 1.2, 1.0, 0.0]
data = [0.0, 12.0, 4.0, 4.0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 0.0]
result = scale_model_sfs_to_data(model, data, verbose=False)
print(round(result["theta"], 6))
print([None if v is None else round(v, 3) for v in result["scaled_spectrum"]])
```

**Example output:**

```text
1.032028
[None, 10.32, 5.16, 3.406, 2.58, 2.064, 1.754, 1.445, 1.238, 1.032, None]
```

## 3) `fit_demographic_model_natural_scale`

- **Purpose:** Fit a demographic model by BFGS on the parameters themselves (not their logarithms).
- **Method:** Calls `Inference.optimize`. Every `method` value of the existing `fit_demographic_model` optimizes log parameters, so no parameter reaches this function. The log-likelihood at the fitted parameters is given by `ll_multinom` (or `ll`).
- **API:** `fit_demographic_model_natural_scale(data[required], *, model[required], p0[required], lower_bound[optional], upper_bound[optional], multinom[optional], grid_points[optional], maxiter[optional], verbose[optional])`
- **Input:** `data` is the observed spectrum of one or two populations (a nested list). `model` names a model of that family, `p0` holds positive starting values, and `lower_bound` and `upper_bound` have one entry per parameter. With `multinom=True` theta is rescaled at every step and is not part of `p0`.
- **Output:** A dict with `params`, `log_likelihood`, `function_calls`, `converged` (BFGS warnflag 0) and the options used.

**Example code:**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_natural_scale

result = fit_demographic_model_natural_scale(
    data,
    model="two_epoch",
    p0=[1.5, 0.7],
    lower_bound=[1e-3, 1e-3],
    upper_bound=[100, 10],
    maxiter=40,
    verbose=False,
)
print([round(v, 2) for v in result["params"]], round(result["log_likelihood"], 2))
```

**Example output:**

```text
[2.0, 0.5] -33.55
```

## 4) `fit_demographic_model_constrained`

- **Purpose:** Fit a demographic model under linear parameter constraints.
- **Method:** Calls `Inference.optimize_cons` (SLSQP). Each constraint is `{"coefficients": [...], "offset": b}` for the value `c . params + b`. An equality constraint must be 0 at the solution and an inequality constraint must be at least 0. For example `nu1 = 2 * nu2` is `{"coefficients": [1, -2]}`. The existing `fit_demographic_model` has no constraint argument.
- **API:** `fit_demographic_model_constrained(data[required], *, model[required], p0[required], lower_bound[optional], upper_bound[optional], equality_constraints[optional], inequality_constraints[optional], multinom[optional], grid_points[optional], maxiter[optional], verbose[optional])`
- **Input:** `data` is the observed spectrum of one or two populations (a nested list). `model` names a model of that family, `p0` holds positive starting values, and `lower_bound` and `upper_bound` have one entry per parameter. With `multinom=True` theta is rescaled at every step and is not part of `p0`. `equality_constraints` and `inequality_constraints` are lists of constraint mappings, and at least one is required.
- **Output:** A dict with `params`, `log_likelihood`, `equality_values` and `inequality_values` (the constraint values at the solution), `iterations` and `converged`.

**Example code:**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_constrained

result = fit_demographic_model_constrained(
    data,
    model="two_epoch",
    p0=[1.5, 0.7],
    lower_bound=[1e-3, 1e-3],
    upper_bound=[100, 10],
    equality_constraints=[{"coefficients": [1, -4]}],
    inequality_constraints=[{"coefficients": [1, -1]}],
    maxiter=40,
    verbose=False,
)
nu, time = result["params"]
print(round(nu, 2), round(time, 2), round(nu - 4 * time, 6), result["converged"])
```

**Example output:**

```text
2.0 0.5 0.0 True
```

## 5) `fit_demographic_model_grid`

- **Purpose:** Fit a demographic model by exhaustive search over a parameter grid.
- **Method:** Calls `Inference.optimize_grid`. `parameter_ranges` gives `[start, stop, num_points]` for each model parameter. The grid holds evenly spaced values from `start` to `stop` inclusive, and every combination (at most 2000) is evaluated. It is useful for finding starting values before a local fit. The existing `fit_demographic_model` has no grid search.
- **API:** `fit_demographic_model_grid(data[required], *, model[required], parameter_ranges[required], multinom[optional], grid_points[optional], verbose[optional])`
- **Input:** `data`, `model`, `multinom` and `grid_points` as above. `parameter_ranges` must have one entry per model parameter, with `start` below `stop` and `num_points` at least 2.
- **Output:** A dict with `params` (the best grid point), `log_likelihood`, `num_evaluations`, `parameter_ranges` and the options used.

**Example code:**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_grid

result = fit_demographic_model_grid(
    data,
    model="two_epoch",
    parameter_ranges=[[1.0, 3.0, 3], [0.2, 1.0, 5]],
    verbose=False,
)
print([round(v, 2) for v in result["params"]], result["num_evaluations"])
print(round(result["log_likelihood"], 2))
```

**Example output:**

```text
[2.0, 0.6] 15
-33.77
```

## 6) `fit_demographic_model_nlopt`

- **Purpose:** Fit a demographic model with a derivative-free NLopt algorithm.
- **Method:** Calls `Inference.opt` (the same object as `NLopt_mod.opt`). `algorithm` is `LN_BOBYQA` (default), `LN_COBYLA`, `LN_NELDERMEAD` or `LN_SBPLX`, and the `nlopt` package is required. The search stops after `maxeval` model evaluations, or when the log-likelihood changes by less than `ftol_abs` or the parameters by less than `xtol_abs`. Upstream's `log_opt` option returns wrong parameters, so it is not exposed.
- **API:** `fit_demographic_model_nlopt(data[required], *, model[required], p0[required], lower_bound[optional], upper_bound[optional], algorithm[optional], multinom[optional], grid_points[optional], maxeval[optional], ftol_abs[optional], xtol_abs[optional], verbose[optional])`
- **Input:** `data` is the observed spectrum of one or two populations (a nested list). `model` names a model of that family, `p0` holds positive starting values, and `lower_bound` and `upper_bound` have one entry per parameter. With `multinom=True` theta is rescaled at every step and is not part of `p0`.
- **Output:** A dict with `params`, `log_likelihood`, `optimizer_value` (the optimum NLopt reports), `algorithm` and the options used.

**Example code:**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_nlopt

result = fit_demographic_model_nlopt(
    data,
    model="two_epoch",
    p0=[1.5, 0.7],
    lower_bound=[1e-3, 1e-3],
    upper_bound=[100, 10],
    algorithm="LN_SBPLX",
    maxeval=300,
    verbose=False,
)
print(result["algorithm"], [round(v, 2) for v in result["params"]], round(result["log_likelihood"], 2))
```

**Example output:**

```text
LN_SBPLX [2.0, 0.5] -33.55
```

## 7) `count_tree_topologies`

- **Purpose:** Count the embedded topologies of named sample sets in one chosen tree.
- **Method:** Calls `tskit.combinatorics.tree_count_topologies` for the single tree given by `tree_index`. The existing `count_sample_topologies` evaluates every tree and has no `tree_index` argument; both give the same counts for the same tree (checked in the tests).
- **API:** `count_tree_topologies(tables[required], sample_sets[required], *, tree_index[optional], max_rows[optional], verbose[optional])`
- **Input:** `tables` is the explicit table mapping of `build_tree_sequence()`. `sample_sets` maps names to lists of sample node IDs, each with at least 2 samples. `tree_index` defaults to 0.
- **Output:** A dict with `tree_index`, `interval` (the tree's `[left, right)`), `sample_sets` and `counts`. The keys of `counts` are set names joined with `+`, and each value is a list of `{shape, label, count}` entries sorted by topology rank (`shape`, `label`).

**Example code:**

```python
TABLES = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 2}, {"time": 3}],
    "edges": [
        {"left": 0, "right": 100, "parent": parent, "child": child}
        for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 20, "ancestral_state": "A"},
        {"position": 30, "ancestral_state": "A"},
    ],
    "mutations": [
        {"site": 0, "node": 0, "derived_state": "G"},
        {"site": 1, "node": 4, "derived_state": "G"},
        {"site": 2, "node": 5, "derived_state": "G"},
    ],
}

from dnakit.population import count_tree_topologies

result = count_tree_topologies(TABLES, {"A": [0, 1], "B": [2, 3]}, verbose=False)
print(result["tree_index"], result["interval"])
for label, entries in result["counts"].items():
    print(label, [entry["count"] for entry in entries])
```

**Example output:**

```text
0 [0.0, 100.0]
A [2]
A+B [4]
B [2]
```

## 8) `export_tree_sequence_nexus`

- **Purpose:** Write a tree sequence as a NEXUS file that includes the TREES block.
- **Method:** Calls `tskit.write_nexus`. With `include_trees=True` every marginal tree is written as a Newick line (named `t<left>^<right>`, taxa `n<node id>`). With `include_alignments=None` the sample haplotypes are added whenever the tree sequence has a site. The existing `export_sample_sequences(file_format="nexus")` fixes `include_trees=False` and writes the alignment only. Trees with several roots cannot be written. The destination must not exist.
- **API:** `export_tree_sequence_nexus(tables[required], destination[required], *, include_trees[optional], include_alignments[optional], precision[optional], missing_data_character[optional], max_rows[optional], verbose[optional])`
- **Input:** `tables` as above. `destination` is a new file path. Options are `include_trees`, `include_alignments` (True, False or None), `precision` (1 to 17 or None) and `missing_data_character` (one character, default `?`).
- **Output:** A dict with `destination`, `bytes`, `num_samples`, `num_trees`, `num_sites`, `has_trees_block` and `has_data_block`.

**Example code:**

```python
TABLES = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 2}, {"time": 3}],
    "edges": [
        {"left": 0, "right": 100, "parent": parent, "child": child}
        for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 20, "ancestral_state": "A"},
        {"position": 30, "ancestral_state": "A"},
    ],
    "mutations": [
        {"site": 0, "node": 0, "derived_state": "G"},
        {"site": 1, "node": 4, "derived_state": "G"},
        {"site": 2, "node": 5, "derived_state": "G"},
    ],
}

import tempfile
from pathlib import Path

from dnakit.population import export_tree_sequence_nexus

with tempfile.TemporaryDirectory() as folder:
    target = Path(folder) / "trees.nex"
    result = export_tree_sequence_nexus(TABLES, target, verbose=False)
    lines = target.read_text().splitlines()
print(result["has_trees_block"], result["has_data_block"], result["num_trees"])
print(*lines[-3:], sep="\n")
```

**Example output:**

```text
True True 1
BEGIN TREES;
  TREE t0^100 = [&R] ((n0:1,n1:1):2,(n2:2,n3:2):1);
END;
```

## 9) `calculate_windowed_divergence`

- **Purpose:** Nucleotide divergence (dxy) per base between two populations in windows of base pairs.
- **Method:** Calls `allel.windowed_divergence`. The existing `calculate_sequence_divergence` gives one value for a whole region and has no window argument.
- **API:** `calculate_windowed_divergence(positions[required], first_counts[required], second_counts[required], *, window_size[required], start[optional], stop[optional], step[optional], is_accessible[optional], verbose[optional])`
- **Input:** `positions` are sorted 1-based positions. Count arrays have one row per variant and one column per allele (at least two). `window_size` is the window length in bases, `start`, `stop` and `step` are optional, and `is_accessible` is a list with one boolean per base from position 1; when given, only accessible bases are counted.
- **Output:** A dict with `window_size`, `step` and `windows`. Each window has `start`, `stop`, `sites`, `dxy` and `callable_bases`. Values that cannot be computed are null.

**Example code:**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_windowed_divergence

result = calculate_windowed_divergence(
    pos, first, second, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], row["callable_bases"], round(row["dxy"], 4))
```

**Example output:**

```text
1 40 4 40 0.0469
41 80 4 40 0.0344
81 120 4 40 0.0344
```

## 10) `calculate_windowed_watterson_theta`

- **Purpose:** Watterson's theta per base in windows of base pairs.
- **Method:** Calls `allel.windowed_watterson_theta`. The existing `calculate_watterson_theta` gives one value for the whole sequence.
- **API:** `calculate_windowed_watterson_theta(positions[required], allele_counts[required], *, window_size[required], start[optional], stop[optional], step[optional], is_accessible[optional], verbose[optional])`
- **Input:** `allele_counts` is the count array of one population; the other arguments are as in `calculate_windowed_divergence`.
- **Output:** A dict with `window_size`, `step` and `windows`. Each window has `start`, `stop`, `sites`, `theta_per_base` and `callable_bases`.

**Example code:**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_windowed_watterson_theta

result = calculate_windowed_watterson_theta(
    pos, first, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], round(row["theta_per_base"], 5))
```

**Example output:**

```text
1 40 4 0.04091
41 80 4 0.02727
81 120 4 0.02727
```

## 11) `calculate_moving_tajimas_d`

- **Purpose:** Tajima's D in moving windows of consecutive variants.
- **Method:** Calls `allel.moving_tajima_d`. The existing `calculate_windowed_tajimas_d` cuts windows by base pairs, which is a different window definition. A window with fewer than `min_sites` variants gives null.
- **API:** `calculate_moving_tajimas_d(allele_counts[required], *, window_size[required], start[optional], stop[optional], step[optional], min_sites[optional], verbose[optional])`
- **Input:** Count arrays with one row per variant. `window_size` is the number of variants per window, `start` and `stop` are variant indices, and `step` defaults to `window_size`.
- **Output:** A dict with `window_size`, `step` and `values` (one per window).

**Example code:**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_moving_tajimas_d

result = calculate_moving_tajimas_d(first, window_size=6, verbose=False)
print(result["window_size"], [round(v, 4) for v in result["values"]])
```

**Example output:**

```text
6 [0.1677, -0.065]
```

## 12) `calculate_weir_cockerham_fst_per_variant`

- **Purpose:** Weir-Cockerham variance components and Fst for every variant.
- **Method:** Calls `allel.weir_cockerham_fst` and sums over alleles. `a` is the among-population component, `b` the among-individuals-within-population component and `c` the within-individual component, and `fst = a / (a + b + c)`. The existing `calculate_fst(method="weir_cockerham")` returns only the jackknifed average and has no per-variant output.
- **API:** `calculate_weir_cockerham_fst_per_variant(genotypes[required], populations[required], *, verbose[optional])`
- **Input:** `genotypes` has shape (variants, samples, ploidy), with -1 for missing calls. `populations` lists the sample indices of each population (at least two, not overlapping).
- **Output:** A dict with `variants`, `populations` (sample counts), the per-variant `a`, `b`, `c` and `fst`, and the overall `fst_ratio_of_sums`.

**Example code:**

```python
genotypes = [
    [[int((v * (s + 1) + h) % 3 == 0) for h in range(2)] for s in range(4)] for v in range(12)
]
groups = [[0, 1], [2, 3]]

from dnakit.variants import calculate_weir_cockerham_fst_per_variant

result = calculate_weir_cockerham_fst_per_variant(genotypes, groups, verbose=False)
print(result["variants"], result["populations"])
print([round(v, 3) for v in result["fst"][:6]])
print(round(result["fst_ratio_of_sums"], 4))
```

**Example output:**

```text
12 [2, 2]
[0.0, -0.333, 0.0, 0.0, -0.333, 0.0]
-0.0909
```

## 13) `calculate_windowed_hudson_fst`

- **Purpose:** Hudson's Fst (ratio of sums) in windows of base pairs.
- **Method:** Calls `allel.windowed_hudson_fst`. The existing `calculate_fst` has no window argument.
- **API:** `calculate_windowed_hudson_fst(positions[required], first_counts[required], second_counts[required], *, window_size[required], start[optional], stop[optional], step[optional], verbose[optional])`
- **Input:** `positions`, the count arrays of two populations and the window arguments, as in `calculate_windowed_divergence` (without `is_accessible`).
- **Output:** A dict with `window_size`, `step` and `windows`. Each window has `start`, `stop`, `sites` and `fst`.

**Example code:**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_windowed_hudson_fst

result = calculate_windowed_hudson_fst(
    pos, first, second, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], round(row["fst"], 4))
```

**Example output:**

```text
1 40 4 0.1111
41 80 4 0.2121
81 120 4 0.0303
```

## 14) `calculate_windowed_weir_cockerham_fst`

- **Purpose:** Weir-Cockerham Fst in windows of base pairs.
- **Method:** Calls `allel.windowed_weir_cockerham_fst`. The existing `calculate_fst` has no window argument.
- **API:** `calculate_windowed_weir_cockerham_fst(positions[required], genotypes[required], populations[required], *, window_size[required], start[optional], stop[optional], step[optional], verbose[optional])`
- **Input:** `positions`, `genotypes` and `populations` (as in the per-variant function), and the window arguments.
- **Output:** A dict with `window_size`, `step` and `windows`. Each window has `start`, `stop`, `sites` and `fst`.

**Example code:**

```python
genotypes = [
    [[int((v * (s + 1) + h) % 3 == 0) for h in range(2)] for s in range(4)] for v in range(12)
]
groups = [[0, 1], [2, 3]]
pos = list(range(10, 130, 10))

from dnakit.variants import calculate_windowed_weir_cockerham_fst

result = calculate_windowed_weir_cockerham_fst(
    pos, genotypes, groups, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], round(row["fst"], 4))
```

**Example output:**

```text
1 40 4 -0.0667
41 80 4 -0.1429
81 120 4 -0.0667
```

## 15) `calculate_moving_hudson_fst`

- **Purpose:** Hudson's Fst in moving windows of consecutive variants.
- **Method:** Calls `allel.moving_hudson_fst`. The existing `calculate_fst` has no window argument.
- **API:** `calculate_moving_hudson_fst(first_counts[required], second_counts[required], *, window_size[required], start[optional], stop[optional], step[optional], verbose[optional])`
- **Input:** Count arrays with one row per variant. `window_size` is the number of variants per window, `start` and `stop` are variant indices, and `step` defaults to `window_size`.
- **Output:** A dict with `window_size`, `step` and `values`.

**Example code:**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_moving_hudson_fst

result = calculate_moving_hudson_fst(first, second, window_size=4, step=2, verbose=False)
print(result["step"], [round(v, 4) for v in result["values"]])
```

**Example output:**

```text
2 [0.1111, 0.0667, 0.2121, 0.1282, 0.0303]
```

## 16) `calculate_moving_weir_cockerham_fst`

- **Purpose:** Weir-Cockerham Fst in moving windows of consecutive variants.
- **Method:** Calls `allel.moving_weir_cockerham_fst`. The existing `calculate_fst` has no window argument.
- **API:** `calculate_moving_weir_cockerham_fst(genotypes[required], populations[required], *, window_size[required], start[optional], stop[optional], step[optional], verbose[optional])`
- **Input:** `genotypes`, `populations` and the moving-window arguments.
- **Output:** A dict with `window_size`, `step` and `values`.

**Example code:**

```python
genotypes = [
    [[int((v * (s + 1) + h) % 3 == 0) for h in range(2)] for s in range(4)] for v in range(12)
]
groups = [[0, 1], [2, 3]]

from dnakit.variants import calculate_moving_weir_cockerham_fst

result = calculate_moving_weir_cockerham_fst(genotypes, groups, window_size=4, step=4, verbose=False)
print(result["step"], [round(v, 4) for v in result["values"]])
```

**Example output:**

```text
4 [-0.0667, -0.1429, -0.0667]
```

## 17) `calculate_blockwise_patterson_d`

- **Purpose:** Patterson's D (ABBA-BABA) with a block-jackknife standard error and z-score.
- **Method:** Calls `allel.blockwise_patterson_d` (`allel.average_patterson_d` is the same object). The four populations are A, B, C and D, and each block holds `block_length` consecutive variants. The existing `calculate_patterson_d` gives only the ratio of sums and no standard error.
- **API:** `calculate_blockwise_patterson_d(first_counts[required], second_counts[required], third_counts[required], fourth_counts[required], *, block_length[required], verbose[optional])`
- **Input:** `first_counts` to `fourth_counts` are the count arrays of the four populations (same number of variants). `block_length` is the block size in variants.
- **Output:** A dict with `d`, `standard_error`, `z_score`, `block_values`, `jackknife_values`, `variants` and `block_length`.

**Example code:**

```python
def counts(m, c):
    return [[4 - ((v * m + c) % 5), (v * m + c) % 5] for v in range(12)]


a, b, c, d = counts(1, 0), counts(2, 1), counts(1, 2), counts(3, 0)

from dnakit.variants import calculate_blockwise_patterson_d

result = calculate_blockwise_patterson_d(a, b, c, d, block_length=4, verbose=False)
print(round(result["d"], 4), round(result["standard_error"], 4), round(result["z_score"], 4))
print([round(v, 4) for v in result["block_values"]])
```

**Example output:**

```text
-0.0442 0.1801 -0.2454
[0.0, 0.3404, -0.4364]
```

## 18) `calculate_moving_patterson_d`

- **Purpose:** Patterson's D in moving windows of consecutive variants.
- **Method:** Calls `allel.moving_patterson_d`. The existing `calculate_patterson_d` has no window argument.
- **API:** `calculate_moving_patterson_d(first_counts[required], second_counts[required], third_counts[required], fourth_counts[required], *, window_size[required], start[optional], stop[optional], step[optional], verbose[optional])`
- **Input:** The count arrays of four populations and the moving-window arguments.
- **Output:** A dict with `window_size`, `step` and `values`.

**Example code:**

```python
def counts(m, c):
    return [[4 - ((v * m + c) % 5), (v * m + c) % 5] for v in range(12)]


a, b, c, d = counts(1, 0), counts(2, 1), counts(1, 2), counts(3, 0)

from dnakit.variants import calculate_moving_patterson_d

result = calculate_moving_patterson_d(a, b, c, d, window_size=4, verbose=False)
print(result["window_size"], [round(v, 4) for v in result["values"]])
```

**Example output:**

```text
4 [0.0, 0.3404, -0.4364]
```

## 19) `calculate_moving_garud_h`

- **Purpose:** Garud's H1, H12, H123 and H2/H1 in moving windows of consecutive variants.
- **Method:** Calls `allel.moving_garud_h`; every window is computed on its own. The existing `calculate_garud_h` computes once for the whole array.
- **API:** `calculate_moving_garud_h(haplotypes[required], *, window_size[required], start[optional], stop[optional], step[optional], verbose[optional])`
- **Input:** `haplotypes` is a 0/1 array of shape (variants, haplotypes) with at least two haplotypes. The window arguments are as above.
- **Output:** A dict with `window_size`, `step` and the lists `h1`, `h12`, `h123` and `h2_h1` (one value per window).

**Example code:**

```python
haplotypes = [[int((v * (h + 1)) % 3 == 0) for h in range(6)] for v in range(12)]

from dnakit.variants import calculate_moving_garud_h

result = calculate_moving_garud_h(haplotypes, window_size=6, verbose=False)
print([round(v, 4) for v in result["h1"]])
print([round(v, 4) for v in result["h12"]])
print([round(v, 4) for v in result["h2_h1"]])
```

**Example output:**

```text
[0.5556, 0.5556]
[1.0, 1.0]
[0.2, 0.2]
```

## 20) `calculate_rogers_huff_r_between`

- **Purpose:** Rogers-Huff r between every variant of one set and every variant of another.
- **Method:** Calls `allel.rogers_huff_r_between`. The existing `calculate_rogers_huff_r` pairs variants inside one set only. Each set may hold at most 500 variants.
- **API:** `calculate_rogers_huff_r_between(alt_counts_a[required], alt_counts_b[required], *, verbose[optional])`
- **Input:** `alt_counts_a` and `alt_counts_b` are alternate-allele dosages of shape (variants, samples) with 0, 1, 2 and -1 for missing. They must cover the same samples, at least two.
- **Output:** A dict with `variants_a`, `variants_b` and the matrix `r` (rows for the first set, columns for the second; r is null for a constant variant).

**Example code:**

```python
from dnakit.variants import calculate_rogers_huff_r_between

result = calculate_rogers_huff_r_between(
    [[0, 1, 2, 1]], [[0, 1, 2, 1], [2, 1, 0, 1], [1, 1, 1, 1]], verbose=False
)
print(result["variants_a"], result["variants_b"])
print(result["r"])
```

**Example output:**

```text
1 3
[[1.0, -1.0, None]]
```

## 21) `calculate_randomized_pca`

- **Purpose:** Randomized principal component analysis of genotype dosages.
- **Method:** Drops monomorphic variants and calls `allel.randomized_pca`. The random state is fixed by `seed`, so the same `seed` gives the same result. The randomized solver approximates the exact PCA of `calculate_pca`, which has no solver argument. `scaler` is `patterson` (default), `standard` or `none`.
- **API:** `calculate_randomized_pca(alt_counts[required], *, n_components[optional], iterated_power[optional], seed[optional], scaler[optional], verbose[optional])`
- **Input:** `alt_counts` has shape (variants, samples) with dosages 0, 1, 2 and no missing values. `n_components` is at most 50 and cannot exceed the number of samples or polymorphic variants. `iterated_power` and `seed` control the solver.
- **Output:** A dict with `coordinates` (one row per sample), `explained_variance_ratio`, `polymorphic_variants`, `n_components`, `seed` and `scaler`.

**Example code:**

```python
dosage = [[(v * (s + 1) + s // 2) % 3 for s in range(6)] for v in range(12)]

from dnakit.variants import calculate_randomized_pca

first = calculate_randomized_pca(dosage, n_components=2, seed=3, verbose=False)
again = calculate_randomized_pca(dosage, n_components=2, seed=3, verbose=False)
print(first["polymorphic_variants"], len(first["coordinates"]), len(first["coordinates"][0]))
print([round(v, 2) for v in first["explained_variance_ratio"]])
print(first == again)
```

**Example output:**

```text
12 6 2
[0.5, 0.3]
True
```

## 22) `calculate_site_frequency_spectrum_from_counts`

- **Purpose:** Site frequency spectrum from biallelic allele counts, optionally folded and scaled.
- **Method:** The unfolded form calls `allel.sfs` on the alternate counts (`sfs_scaled` with `scale=True`); `fold=True` calls `sfs_folded` (`sfs_folded_scaled` with `scale=True`). `spectrum[i]` counts variants with `i` derived alleles, and a folded spectrum has `sample_size // 2 + 1` classes. The existing `calculate_site_frequency_spectrum` reads variant records and has no fold or scale option.
- **API:** `calculate_site_frequency_spectrum_from_counts(allele_counts[required], *, fold[optional], scale[optional], sample_size[optional], verbose[optional])`
- **Input:** `allele_counts` has shape (variants, 2): reference and alternate (derived) counts. `sample_size` defaults to the largest total over variants and may not be below it.
- **Output:** A dict with `variants`, `sample_size`, `folded`, `scaled` and `spectrum` (floats when scaled, integers otherwise).

**Example code:**

```python
from dnakit.variants import calculate_site_frequency_spectrum_from_counts

counts = [[3, 1], [2, 2], [0, 4], [4, 0], [1, 3]]
for fold, scale in [(False, False), (True, False), (False, True), (True, True)]:
    result = calculate_site_frequency_spectrum_from_counts(
        counts, fold=fold, scale=scale, verbose=False
    )
    print(fold, scale, result["sample_size"], result["spectrum"])
```

**Example output:**

```text
False False 4 [1, 1, 1, 1, 1]
True False 4 [2, 2, 1]
False True 4 [0.0, 1.0, 2.0, 3.0, 4.0]
True True 4 [0.0, 1.5, 1.0]
```

## 23) `calculate_joint_site_frequency_spectrum_from_counts`

- **Purpose:** Joint site frequency spectrum of two populations from biallelic allele counts.
- **Method:** The unfolded form calls `allel.joint_sfs` (`joint_sfs_scaled` with `scale=True`); `fold=True` calls `joint_sfs_folded` (`joint_sfs_folded_scaled` with `scale=True`). Each population is then folded to its own minor allele, which differs from `fold_joint_site_frequency_spectrum`, which folds the pooled sample. The existing one-population spectrum wrappers cannot return a joint spectrum.
- **API:** `calculate_joint_site_frequency_spectrum_from_counts(first_counts[required], second_counts[required], *, fold[optional], scale[optional], first_sample_size[optional], second_sample_size[optional], verbose[optional])`
- **Input:** `first_counts` and `second_counts` are (variants, 2) counts of the same variants in two populations. `first_sample_size` and `second_sample_size` default to the largest total of each.
- **Output:** A dict with `variants`, both sample sizes, `folded`, `scaled` and the two-dimensional `spectrum`. `spectrum[i][j]` corresponds to `i` derived alleles in the first and `j` in the second population.

**Example code:**

```python
from dnakit.variants import calculate_joint_site_frequency_spectrum_from_counts

first = [[3, 1], [2, 2], [0, 4]]
second = [[4, 0], [1, 3], [2, 2]]
plain = calculate_joint_site_frequency_spectrum_from_counts(first, second, verbose=False)
folded = calculate_joint_site_frequency_spectrum_from_counts(first, second, fold=True, verbose=False)
print(plain["first_sample_size"], plain["second_sample_size"])
for row in plain["spectrum"]:
    print(row)
print(folded["spectrum"])
```

**Example output:**

```text
4 4
[0, 0, 0, 0, 0]
[1, 0, 0, 0, 0]
[0, 0, 0, 1, 0]
[0, 0, 0, 0, 0]
[0, 0, 1, 0, 0]
[[0, 0, 1], [1, 0, 0], [0, 1, 0]]
```

## 24) `phase_trio_step_by_transmission`

- **Purpose:** One step of trio phasing by Mendelian transmission.
- **Method:** `step="progeny"` calls `allel.phase_progeny_by_transmission` and phases the child where the parents allow it. `step="parents"` calls `allel.phase_parents_by_transmission` and phases the parents from the phased child within windows of `window_size` variants; its input must be the `genotypes` and `is_phased` returned by the first step. The existing `phase_trio_by_transmission` always runs both steps and cannot take phased input.
- **API:** `phase_trio_step_by_transmission(genotypes[required], *, step[optional], window_size[optional], is_phased[optional], verbose[optional])`
- **Input:** `genotypes` has shape (variants, 3, 2) in the order mother, father, child, with -1 for missing. `step` is `progeny` or `parents`; the `parents` step also needs `window_size` and a boolean array `is_phased` of shape (variants, 3).
- **Output:** A dict with `step`, `variants`, `window_size`, the phased `genotypes` and `is_phased`.

**Example code:**

```python
trio = [
    [[(v + s) % 2, (v * (s + 1) + 1) % 2] for s in range(3)] for v in range(12)
]

from dnakit.variants import phase_trio_step_by_transmission

one = phase_trio_step_by_transmission(trio, verbose=False)
two = phase_trio_step_by_transmission(
    one["genotypes"], step="parents", window_size=4, is_phased=one["is_phased"], verbose=False
)
print(sum(row[2] for row in one["is_phased"]), sum(sum(row) for row in one["is_phased"]))
print(sum(sum(row) for row in two["is_phased"]), two["window_size"])
```

**Example output:**

```text
6 6
18 4
```

## 25) `read_vcf_field_table`

- **Purpose:** Read the variant fields of a local VCF into a flat table.
- **Method:** Calls `allel.vcf_to_dataframe` (`allel.vcf_to_recarray` is the NumPy record-array form of the same table). `ALT` and multi-value INFO fields are spread over `alt_number` columns (`ALT_1`, `ALT_2`, ...). `fields` defaults to the standard fixed columns and may name fields such as `variants/DP`. The existing `read_vcf` returns nested records with genotypes and has no field selection or `alt_number`.
- **API:** `read_vcf_field_table(path[required], *, fields[optional], alt_number[optional], max_rows[optional], verbose[optional])`
- **Input:** `path` is a local VCF. `fields` is a list of allel field names or None, `alt_number` is between 1 and 20, and `max_rows` caps the records (at most 100000).
- **Output:** A dict with `columns`, `num_rows` and `rows` (one dict per record, missing values as null).

**Example file:** `x.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	10	rs1	A	G	50	PASS	DP=20;AF=0.5	GT	0/1	1/1
chr1	20	.	C	T,A	.	PASS	DP=7;AF=0.2,0.1	GT	0|1	0/0
chr1	30	rs3	G	C	99	PASS	DP=12	GT	0/0	0/1
```

**Example code:**

```python
from dnakit.variants import read_vcf_field_table

result = read_vcf_field_table(
    "x.vcf", fields=["CHROM", "POS", "REF", "ALT", "variants/DP"], alt_number=2, verbose=False
)
print(result["columns"], result["num_rows"])
for row in result["rows"]:
    print(row)
```

**Example output:**

```text
['CHROM', 'POS', 'REF', 'ALT_1', 'ALT_2', 'DP'] 3
{'CHROM': 'chr1', 'POS': 10, 'REF': 'A', 'ALT_1': 'G', 'ALT_2': None, 'DP': 20}
{'CHROM': 'chr1', 'POS': 20, 'REF': 'C', 'ALT_1': 'T', 'ALT_2': 'A', 'DP': 7}
{'CHROM': 'chr1', 'POS': 30, 'REF': 'G', 'ALT_1': 'C', 'ALT_2': None, 'DP': 12}
```

## 26) `vcf_to_csv_table`

- **Purpose:** Write the variant fields of a local VCF to a new CSV file.
- **Method:** Calls `allel.vcf_to_csv`; the columns are those of `read_vcf_field_table`. The output must end in `.csv` and must not exist, and a partly written file is removed on error. The existing `vcf_to_npz_archive` and `vcf_to_hdf5_archive` write arrays, not a table file.
- **API:** `vcf_to_csv_table(source[required], output[required], *, fields[optional], alt_number[optional], verbose[optional])`
- **Input:** `source` is a local VCF, `output` is a new CSV path, and `fields` and `alt_number` are as in `read_vcf_field_table`.
- **Output:** A dict with `source`, `output`, `columns`, `num_rows` and `bytes`.

**Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_to_csv_table

with tempfile.TemporaryDirectory() as folder:
    target = Path(folder) / "table.csv"
    result = vcf_to_csv_table("x.vcf", str(target), fields=["CHROM", "POS", "REF", "ALT"], verbose=False)
    lines = target.read_text().splitlines()
print(result["columns"], result["num_rows"])
print(*lines, sep="\n")
```

**Example output:**

```text
['CHROM', 'POS', 'REF', 'ALT_1', 'ALT_2', 'ALT_3'] 3
CHROM,POS,REF,ALT_1,ALT_2,ALT_3
chr1,10,A,G,,
chr1,20,C,T,A,
chr1,30,G,C,,
```
