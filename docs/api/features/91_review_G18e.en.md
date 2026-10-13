# Review wrappers for dadi and tskit (group G18e, doc 91)

This page records the group G18e review of the functions marked `other` in dadi and tskit. The dadi functions that are DNA-level and can run offline are wrapped in `dnakit.population`; the rest keep a precise status and reason. None of the seven tskit functions is a DNA-level operation, so `dnakit.population` exports nothing.

The wrappers are exported from `dnakit.population`, and the examples import them from that package. Dependencies are optional and are imported only when a function is called. Output files are never overwritten; model functions are chosen by name from fixed lists and no arbitrary code runs; the global random state is restored after each call.

## Example data

The files below are written into the working directory before the examples run. The other example data are generated in the code from a dadi model with Poisson noise and a fixed seed.

**Example file：** `variants.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3	s4
chr1	100	.	A	G	.	PASS	AA=A	GT	0|0	0|1	1|1	0|0
chr1	200	.	C	T	.	PASS	AA=C	GT	0|1	1|1	0|0	0|0
chr1	300	.	G	A	.	PASS	AA=G	GT	0|0	0|0	0|1	1|1
chr1	400	.	T	C	.	PASS	AA=T	GT	1|0	0|1	0|0	1|1
```

**Example file：** `popinfo.txt`

```text
s1 popA
s2 popA
s3 popB
s4 popB
```

**Example file：** `slim_p1.txt`

```text
#OUT: 1 1 A
Mutations:
0 0 m1 100 0.1
1 2 m1 250 0.2
2 5 m2 400 0.0
Genomes:
p1:0 A 0 1
p1:1 A 1
p1:2 A
```

**Example file：** `slim_p2.txt`

```text
#OUT: 1 1 A
Mutations:
0 0 m1 100 0.1
1 3 m1 300 0.2
Genomes:
p2:0 A 0
p2:1 A 1
p2:2 A 0 1
```

## 1) `dnakit.population.total_substitution_rate()`

- **Purpose:** Total instantaneous substitution rate of a DNA Markov model.
- **Method:** Calls dadi's `Misc.total_instantaneous_rate()`. The diagonal is ignored; the total is the sum over states of frequency times exit rate.
- **API:** `dnakit.population.total_substitution_rate(rate_matrix, stationary, *, verbose=True)`
- **Input:** `rate_matrix` is square over 4, 16 or 64 states (bases, dinucleotides or trinucleotides in dadi's CGTA order), with `rate_matrix[i][j]` the rate from state i to j; `stationary` gives each state's frequency and must sum to 1.
- **Output:** A dictionary with `num_states`, `total_rate` and `provenance`.
- **Example code:**

```python
from dnakit.population import total_substitution_rate

jc = [[0.0 if i == j else 1 / 3 for j in range(4)] for i in range(4)]
result = total_substitution_rate(jc, [0.25] * 4)
print(result["num_states"], round(result["total_rate"], 6))
```

- **Example output:**

```text
4 1.0
```

## 2) `dnakit.population.ms_command_line()`

- **Purpose:** Build an `ms` command line for coalescent sequence simulation; the text is built only, and `ms` is not run.
- **Method:** Calls dadi's `Misc.ms_command()`. When recombination is set, the recombination sites default to `10 * theta`.
- **API:** `dnakit.population.ms_command_line(theta, sample_sizes, core, iterations, *, recombination=0.0, recombination_sites=None, seeds=None, verbose=True)`
- **Input:** `theta` is positive; `sample_sizes` gives one haploid sample size per population; `core` is a one-line demography fragment; `iterations` is the number of simulations (at most 100000); `seeds` is three integers or None.
- **Output:** A dictionary with `command`, `theta`, `sample_sizes`, `iterations` and `provenance`.
- **Example code:**

```python
from dnakit.population import ms_command_line

result = ms_command_line(2.0, [3, 4], "-eN 0.5 2.0", 10, recombination=1.0, seeds=(1, 2, 3))
print(result["command"])
```

- **Example output:**

```text
ms 7 10 -t 2.000000 -I 2 3 4 -eN 0.5 2.0 -r 1.000000 20 -seeds 1 2 3
```

## 3) `dnakit.population.ms_core_command_2d()`

- **Purpose:** Build the `ms` core command (the demography part) of a two-population model.
- **Method:** Calls dadi's `Demographics2D.*_mscore()`. Times and rates use dadi's scaling, not coalescent units.
- **API:** `dnakit.population.ms_core_command_2d(model, params, *, verbose=True)`
- **Input:** `model` is `split_mig` (parameters nu1, nu2, T, m), `IM` (s, nu1, nu2, T, m12, m21) or `IM_pre` (nuPre, TPre, s, nu1, nu2, T, m12, m21); `s` must lie strictly between 0 and 1.
- **Output:** A dictionary with `model`, `core` and `provenance`.
- **Example code:**

```python
from dnakit.population import ms_core_command_2d

print(ms_core_command_2d("split_mig", [1.0, 2.0, 0.5, 0.25])["core"])
```

- **Example output:**

```text
-n 1 1.000000 -n 2 2.000000 -ma x 0.500000 0.500000 x -ej 0.250000 2 1 -en 0.250000 1 1
```

## 4) `dnakit.population.expected_neutral_sfs_2d()`

- **Purpose:** Expected neutral spectrum for two populations that never diverge.
- **Method:** Calls dadi's `Demographics2D.snm()`. The populations never diverge, so each is a constant-size neutral sample.
- **API:** `dnakit.population.expected_neutral_sfs_2d(sample_sizes, pts, *, verbose=True)`
- **Input:** `sample_sizes` gives two sample sizes (each at least 2, at most 200); `pts` is the number of integration grid points (20 to 400).
- **Output:** A dictionary with `sample_sizes`, `pts`, `spectrum` and `provenance`; corners that dadi masks are None.
- **Example code:**

```python
from dnakit.population import expected_neutral_sfs_2d

result = expected_neutral_sfs_2d((3, 3), 20)
print(result["sample_sizes"], result["spectrum"][0][0], round(result["spectrum"][0][1], 6))
```

- **Example output:**

```text
[3, 3] None 0.552742
```

## 5) `dnakit.population.trinucleotide_misid_table()`

- **Purpose:** Write the 1 - fux table used to correct ancestral misidentification.
- **Method:** Calls dadi's `Misc.make_fux_table()`, which uses a matrix exponential. The input rate matrix is copied, so the caller's array is not changed. The output file must be new.
- **API:** `dnakit.population.trinucleotide_misid_table(path, *, expected_substitutions, rate_matrix, trinucleotide_frequencies, verbose=True)`
- **Input:** `path` is a new file; `expected_substitutions` is positive and is the expected number of substitutions per site between ingroup and outgroup; `rate_matrix` is the 64 x 64 trinucleotide rate matrix in dadi's CGTA order; `trinucleotide_frequencies` maps the 64 trinucleotides over A, C, G and T to frequencies.
- **Output:** A dictionary with `path`, `rows` (always 256), `expected_substitutions` and `provenance`.
- **Example code:**

```python
from itertools import product
from pathlib import Path
from dnakit.population import trinucleotide_misid_table

keys = ["".join(p) for p in product("ACGT", repeat=3)]
rates = [[0.0 if i == j else 1 / 63 for j in range(64)] for i in range(64)]
Path("misid_table.txt").unlink(missing_ok=True)
result = trinucleotide_misid_table("misid_table.txt", expected_substitutions=0.1, rate_matrix=rates,
    trinucleotide_frequencies={key: 1 / 64 for key in keys})
print(result["rows"], Path(result["path"]).read_text().splitlines()[1])
```

- **Example output:**

```text
256 CCC G 0.001626
```

## 6) `dnakit.population.slim_output_data_dict()`

- **Purpose:** Read SLiM sample outputs, one file per population, into a dadi SNP dictionary.
- **Method:** Calls dadi's `Misc.dd_from_SLiM_files()`. Populations are numbered 0, 1, ... in file order; each site is keyed `<prefix><position>.<global mutation id>`; results are sorted by key.
- **API:** `dnakit.population.slim_output_data_dict(paths, *, mutation_types=None, prefix="SLIM_", verbose=True)`
- **Input:** `paths` lists one SLiM output file per population (at most 10); `mutation_types` keeps only the listed SLiM mutation types, and None keeps all.
- **Output:** A dictionary with `populations`, `sample_sizes`, `num_sites` and `sites`; each site's `calls` gives [ancestral count, derived count] per population.
- **Example code:**

```python
from dnakit.population import slim_output_data_dict

result = slim_output_data_dict(["slim_p1.txt", "slim_p2.txt"])
print(result["sample_sizes"], result["num_sites"])
print(result["sites"][0])
print(slim_output_data_dict(["slim_p1.txt", "slim_p2.txt"], mutation_types=["m2"])["num_sites"])
```

- **Example output:**

```text
[3, 3] 4
{'key': 'SLiM_100.0', 'outgroup_allele': 0, 'calls': {'0': [2, 1], '1': [1, 2]}}
1
```

## 7) `dnakit.population.vcf_subsample_bootstrap_spectra()`

- **Purpose:** Bootstrap spectra from a VCF, drawing a random subset of individuals and chromosome chunks each time.
- **Method:** Calls dadi's `Misc.bootstraps_subsample_vcf()`. `seed` fixes both the individual subsets and the chunk draws; the global Python and NumPy random states are saved and restored.
- **API:** `dnakit.population.vcf_subsample_bootstrap_spectra(vcf_path, popinfo_path, *, subsample, pop_ids, num_boot, chunk_size, seed=0, filter_sites=True, polarized=True, verbose=True)`
- **Input:** `subsample` gives the number of diploid individuals kept per population; `pop_ids` lists one to three populations; `num_boot` is at most 200; `chunk_size` is the chunk length in base pairs.
- **Output:** A dictionary with `num_boot`, `seed`, `spectra` (each bootstrap spectrum is a nested list with corners None) and `provenance`.
- **Example code:**

```python
from dnakit.population import vcf_subsample_bootstrap_spectra

result = vcf_subsample_bootstrap_spectra("variants.vcf", "popinfo.txt",
    subsample={"popA": 2, "popB": 2}, pop_ids=["popA", "popB"], num_boot=2, chunk_size=200, seed=3)
print(result["num_boot"], len(result["spectra"]), result["spectra"][0][0])
```

- **Example output:**

```text
2 2 [None, 0.0, 0.0, 1.0, 0.0]
```

## 8) `dnakit.population.parameter_uncertainty_gim()`

- **Purpose:** Parameter standard deviations from the Godambe information matrix (GIM).
- **Method:** Calls dadi's `Godambe.GIM_uncert()`. The model is chosen by name from a fixed list and extrapolated over `grid_points`. With `multinom=True`, theta is estimated too and its uncertainty is the last entry.
- **API:** `dnakit.population.parameter_uncertainty_gim(data, bootstraps, *, model, p0, grid_points=(20, 30, 40), log=False, multinom=True, eps=0.01, verbose=True)`
- **Input:** `data` is a spectrum (nested list, corners may be None); `bootstraps` is a list of spectra with the shape of `data`; `p0` holds positive starting parameters.
- **Output:** A dictionary with `params`, `uncertainties`, `num_boot` and `provenance`.
- **Example code:**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
mean = np.ma.filled(model((1.5, 0.5), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(7)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import parameter_uncertainty_gim

result = parameter_uncertainty_gim(data, boots, model="two_epoch", p0=[1.5, 0.5])
print([round(value, 3) for value in result["uncertainties"]])
```

- **Example output:**

```text
[0.112, 0.109, 6.849]
```

## 9) `dnakit.population.godambe_matrices()`

- **Purpose:** Godambe information, Hessian and variability matrices.
- **Method:** Calls dadi's `Godambe.get_godambe()`, which uses dadi's Poisson log-likelihood. The model is evaluated as given, so include theta among the parameters if it should be estimated. The Godambe information is H J^-1 H.
- **API:** `dnakit.population.godambe_matrices(data, bootstraps, *, model, p0, grid_points=(20, 30, 40), log=False, eps=0.01, verbose=True)`
- **Input:** The same as `parameter_uncertainty_gim()`.
- **Output:** A dictionary with `godambe_information`, `hessian`, `variability` and `provenance`.
- **Example code:**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
mean = np.ma.filled(model((1.5, 0.5), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(7)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import godambe_matrices

result = godambe_matrices(data, boots, model="two_epoch", p0=[1.5, 0.5])
print(len(result["godambe_information"]), round(result["godambe_information"][0][0], 3))
```

- **Example output:**

```text
2 87951.897
```

## 10) `dnakit.population.nested_model_statistic()`

- **Purpose:** Nested-model tests with bootstrap adjustment: the likelihood-ratio (LRT) adjustment factor, the adjusted Wald statistic and the adjusted score statistic.
- **Method:** Calls dadi's `Godambe.LRT_adjust()`, `Wald_stat()` and `score_stat()`. `nested_indices` gives the positions of the nested parameters in the complex model's list. A zero nested value makes the matrices singular, so use a non-zero value.
- **API:** `dnakit.population.nested_model_statistic(data, bootstraps, *, test, model, p0, nested_indices, grid_points=(20, 30, 40), full_params=None, multinom=True, eps=0.01, adjusted_and_original=False, verbose=True)`
- **Input:** `test` is `lrt`, `wald` or `score`; `p0` holds the complex model's parameters with the nested ones at their nested values; the Wald test also needs `full_params`, the complete list of complex-model parameter values.
- **Output:** `lrt` returns `adjustment`; `wald` and `score` return `statistic`, and `adjusted_and_original=True` adds `unadjusted`.
- **Example code:**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.three_epoch)
mean = np.ma.filled(model((2.0, 0.5, 0.2, 0.3), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(3)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import nested_model_statistic

result = nested_model_statistic(data, boots, test="score", model="three_epoch",
    p0=[1.0, 0.5, 0.2, 0.3], nested_indices=[0], adjusted_and_original=True)
print(round(result["statistic"], 6), round(result["unadjusted"], 6))
```

- **Example output:**

```text
0.906362 2.542513
```

## 11) `dnakit.population.fit_demographic_model_to_residual()`

- **Purpose:** Fit a demographic model so that its residual spectrum matches a given target residual spectrum.
- **Method:** Calls dadi's `Inference.optimize_log_resid()`, a BFGS optimizer on log parameters. The number of iterations is capped at `max_iterations`.
- **API:** `dnakit.population.fit_demographic_model_to_residual(data, target_residual, *, model, p0, grid_points=(20, 30, 40), lower_bound=None, upper_bound=None, multinom=True, max_iterations=50, verbose=True)`
- **Input:** `target_residual` is the target residual spectrum, with the same shape as `data` (for example from a synonymous-site fit); `p0` holds positive values.
- **Output:** A dictionary with `params`, `objective`, `converged`, `function_calls`, `gradient_calls` and `provenance`.
- **Example code:**

```python
import numpy as np
import dadi
model = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
mean = np.ma.filled(model((1.5, 0.5), (10,), [20, 30, 40]) * 200, 0.0)
rng = np.random.default_rng(7)
data = rng.poisson(mean).astype(float).tolist()
boots = [rng.poisson(mean).astype(float).tolist() for _ in range(4)]
from dnakit.population import fit_demographic_model_to_residual

target = [value * 0.5 for value in data]
result = fit_demographic_model_to_residual(data, target, model="two_epoch", p0=[0.8, 0.3], max_iterations=2)
print([round(value, 1) for value in result["params"]], result["converged"])
```

- **Example output:**

```text
[4060.9, 29915.2] False
```

## 12) tskit: no wrapped functions

None of the seven tskit functions marked `other` is a DNA-level operation, so this module exports nothing:

| Function | Status | Reason |
|---|---|---|
| `tskit.naturalsize` | not-dna | display helper that formats byte counts |
| `tskit.format_number` | not-dna | display helper that formats numbers |
| `tskit.set_print_options` | not-dna | changes global print options; display only |
| `tskit.jit.numba.jitwrap` | internal | converts a tree sequence for numba; not a DNA-level operation |
| `tskit.util.naturalsize`, `format_number`, `set_print_options` | not-dna | the same functions as above (duplicates) |

## dadi functions not wrapped

- **TwoLocus** and **Triallele**: they do not run in this environment. Triallele uses `np.float`, which NumPy 2 removed (status `broken-upstream`). TwoLocus fails with `dadi.tridiag` missing (dadi 2.4.4 names the module `tridiag_cython`), which is not a NumPy 2 failure.
- **DFE mixture models** need precomputed Cache1D/Cache2D selection caches (`needs-data`); **Demes.Inference.optimize** needs a dadi inference-options file (`needs-data`); **LowPass.compute_cov_dist** needs per-site allelic depth in a data dictionary (`needs-data`).
- The other **LowPass** functions are intermediate steps of the low-coverage correction; the **Demes** graph helpers and output depend on module-level global state; **Numerics.make_anc_state_misid_func** takes a Python callable; **Inference.add_misid_param** is deprecated and fails when called.
