# dadi spectrum operations, statistics, low coverage and diffusion

Functions of the `dnakit.population` area that wrap dadi's `Spectrum` methods (masking, projection, population selection and reordering, log, unfolding, resampling, summary statistics, ms and sfs_code parsing), `Numerics` grid extrapolation, `Misc.fragment_data_dict`, `Godambe` derivatives, the `Integration` of the diffusion density for one to five populations, and the GATK multi-sample low-coverage model of `dadi.LowPass`. Spectra are nested lists with `None` for masked entries. Random operations take a `seed` and restore NumPy's global random state. Requires `dnakit[bio-dadi]` (dadi 2.4.4 was tested). `Inference.optimize_lbfgsb` is not offered: it raises `TypeError` (`fmin_l_bfgs_b` no longer accepts `iprint`) for every input.

## Example data

**Example file:** `toy.ms`

```text
ms 4 2 -t 2.0 -I 2 2 2
123 456 789

//
segsites: 2
positions: 0.1 0.5
01
10
01
10

//
segsites: 1
positions: 0.3
1
0
0
1
```

**Example file:** `snps.txt`

```text
Ingroup	Outgroup	Allele1	A	B	Allele2	A	B	Gene	Position
TCT	TCT	C	4	3	T	2	2	g1	chr1_100
ACA	ACA	C	1	6	A	3	1	g1	chr1_200
GCG	GCG	C	1	6	G	3	1	g2	chr2_50
```

## 1) `dnakit.population.calculate_sfs_statistic()`

- **Purpose:** Population-genetic statistic of a site frequency spectrum.
- **Method:** Calls the `Spectrum` method named by `statistic`: `S` (segregating sites), `Watterson_theta`, `theta_L` (Zeng et al. 2006), `Zengs_E`, `pi` (pairwise differences with the n / (n - 1) correction), `Tajima_D` or `Fst` (Weir and Cockerham, joint spectrum). All but `S` and `Fst` need a one-population spectrum. One function covers seven work-list rows.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.calculate_sfs_statistic(spectrum, statistic, *, folded=False, mask_corners=True)`.
- **Input:** `spectrum` (nested list, None for masked), `statistic`, `folded`, `mask_corners`.
- **Output:** `value` and `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import calculate_sfs_statistic

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
for name in ("S", "Watterson_theta", "pi", "Tajima_D"):
    print(name, round(calculate_sfs_statistic(spectrum, name, verbose=False)["value"], 4))
```

- **Example output:**

```text
S 28.29
Watterson_theta 10.0001
pi 10.0002
Tajima_D 0.0
```

## 2) `dnakit.population.mask_spectrum_corners()`

- **Purpose:** Mask the entries seen in none or in all samples.
- **Method:** Calls `Spectrum.mask_corners`. `Spectrum.unmask_all` is not offered because it raises `IndexError` for every input with NumPy 2.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.mask_spectrum_corners(spectrum, *, folded=False)`.
- **Input:** `spectrum` and `folded`.
- **Output:** `spectrum` with the two corner entries None, `sample_sizes`, `folded`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import mask_spectrum_corners

joint = [[(i + 1) * (j + 2) % 5 + 0.5 for j in range(5)] for i in range(5)]
result = mask_spectrum_corners(joint, verbose=False)
print(result["spectrum"][0][0], result["spectrum"][4][4], result["spectrum"][0][1])
```

- **Example output:**

```text
None None 3.5
```

## 3) `dnakit.population.project_spectrum()`

- **Purpose:** Project a spectrum to smaller sample sizes.
- **Method:** Calls `Spectrum.project`, which resamples hypergeometrically. No sample size may exceed the current one.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.project_spectrum(spectrum, sample_sizes, *, folded=False, mask_corners=True)`.
- **Input:** `spectrum`, `sample_sizes` (one per population), `folded`, `mask_corners`.
- **Output:** The projected `spectrum`, `projected_to`, `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import project_spectrum

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
result = project_spectrum(spectrum, [6], verbose=False)
print(result["sample_sizes"], [None if v is None else round(v, 3) for v in result["spectrum"]])
```

- **Example output:**

```text
[6] [None, 10.0, 5.0, 3.333, 2.5, 2.0, None]
```

## 4) `dnakit.population.subset_spectrum_populations()`

- **Purpose:** Keep some populations of a joint spectrum and sum over the others.
- **Method:** `keep` calls `Spectrum.filter_pops`; `marginalize_over` calls `Spectrum.marginalize` (dadi counts its axes from 0, DNAKit numbers populations from 1 in both). Give exactly one. One function covers two work-list rows.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.subset_spectrum_populations(spectrum, *, keep=None, marginalize_over=None, folded=False, mask_corners=True)`.
- **Input:** A joint `spectrum` and `keep` or `marginalize_over` (population numbers from 1), `folded`, `mask_corners`.
- **Output:** `spectrum`, `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import subset_spectrum_populations

joint = [[(i + 1) * (j + 2) % 5 + 0.5 for j in range(5)] for i in range(5)]
kept = subset_spectrum_populations(joint, keep=[2], verbose=False)
summed = subset_spectrum_populations(joint, marginalize_over=[1], verbose=False)
print(kept["sample_sizes"], kept["spectrum"] == summed["spectrum"])
print([None if v is None else round(v, 2) for v in kept["spectrum"]])
```

- **Example output:**

```text
[4] True
[None, 12.5, 12.5, 2.5, None]
```

## 5) `dnakit.population.log_spectrum()`

- **Purpose:** Natural logarithm of a spectrum.
- **Method:** Calls `Spectrum.log`; entries that are not positive become masked.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.log_spectrum(spectrum, *, folded=False, mask_corners=True)`.
- **Input:** `spectrum`, `folded`, `mask_corners`.
- **Output:** `spectrum` of logarithms, `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import log_spectrum

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
result = log_spectrum(spectrum, verbose=False)
print([None if v is None else round(v, 3) for v in result["spectrum"][:4]])
```

- **Example output:**

```text
[None, 2.303, 1.609, 1.204]
```

## 6) `dnakit.population.reorder_spectrum_populations()`

- **Purpose:** Reorder the populations of a joint spectrum.
- **Method:** Calls `Spectrum.reorder_pops`; `new_order` lists the original population numbers (from 1) in the new order.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.reorder_spectrum_populations(spectrum, new_order, *, folded=False, mask_corners=True)`.
- **Input:** `spectrum`, `new_order`, `folded`, `mask_corners`.
- **Output:** The transposed `spectrum`, `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import reorder_spectrum_populations

joint = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
result = reorder_spectrum_populations(joint, [2, 1], verbose=False)
print(result["sample_sizes"], result["spectrum"][1])
```

- **Example output:**

```text
[2, 3] [2.0, 5.0, 8.0, 11.0]
```

## 7) `dnakit.population.unfold_spectrum()`

- **Purpose:** Unfold a folded spectrum.
- **Method:** Calls `Spectrum.unfold`, which assumes each allele is equally likely to be ancestral: the result is `(spectrum + reversed spectrum) / 2`.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.unfold_spectrum(spectrum, *, mask_corners=True)`.
- **Input:** A folded `spectrum` (minor-allele counts, None for masked), `mask_corners`.
- **Output:** The unfolded `spectrum`, `folded` (False), `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import unfold_spectrum

folded = [None, 4.0, 2.0, 1.0, None, None, None]
result = unfold_spectrum(folded, verbose=False)
print(result["folded"], result["spectrum"])
```

- **Example output:**

```text
False [None, 2.0, 1.0, 1.0, 1.0, 2.0, None]
```

## 8) `dnakit.population.resample_spectrum()`

- **Purpose:** Resample a spectrum by Poisson draws or to a fixed number of sites.
- **Method:** `mode='poisson'` calls `Spectrum.sample` (each entry becomes a Poisson draw with the entry as mean); `mode='fixed_size'` calls `Spectrum.fixed_size_sample` (a multinomial draw of `num_sites` sites). `seed` makes the draw reproducible and NumPy's global random state is restored afterwards. One function covers two work-list rows.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.resample_spectrum(spectrum, mode='poisson', *, num_sites=None, only_nonmasked=False, seed=None, folded=False, mask_corners=True)`.
- **Input:** `spectrum`, `mode`, `num_sites`, `only_nonmasked`, `seed`, `folded`, `mask_corners`.
- **Output:** The resampled `spectrum`, `mode`, `seed`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import resample_spectrum

spectrum = [None] + [round(10 / i, 3) for i in range(1, 10)] + [None]
fixed = resample_spectrum(spectrum, "fixed_size", num_sites=100, seed=3, only_nonmasked=True, verbose=False)
print(sum(v for v in fixed["spectrum"] if v is not None))
again = resample_spectrum(spectrum, "fixed_size", num_sites=100, seed=3, only_nonmasked=True, verbose=False)
print(fixed["spectrum"] == again["spectrum"])
```

- **Example output:**

```text
100.0
True
```

## 9) `dnakit.population.spectrum_from_simulation_file()`

- **Purpose:** Build a spectrum from ms or sfs_code simulation output.
- **Method:** `file_format='ms'` calls `Spectrum.from_ms_file` (populations from the ms command line or `pop_assignments`); `file_format='sfs_code'` calls `Spectrum.from_sfscode_file` with `sites` set to `all`, `syn` or `nonsyn`. `average` averages the replicate runs, otherwise they are summed. One function covers two work-list rows.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.spectrum_from_simulation_file(path, *, file_format='ms', average=True, mask_corners=True, pop_assignments=None, sites='all')`.
- **Input:** `path`, `file_format`, `average`, `mask_corners`, `pop_assignments`, `sites`.
- **Output:** `spectrum`, `sample_sizes`, the ms `command` and `seeds`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import spectrum_from_simulation_file

result = spectrum_from_simulation_file("toy.ms", verbose=False)
print(result["sample_sizes"], result["command"])
print(result["spectrum"])
```

- **Example output:**

```text
[2, 2] ms 4 2 -t 2.0 -I 2 2 2
[[None, 0.0, 0.0], [0.0, 1.5, 0.0], [0.0, 0.0, None]]
```

## 10) `dnakit.population.fragment_snp_data()`

- **Purpose:** Split a dadi SNP file into genomic chunks for block bootstrapping.
- **Method:** Reads the SNP file with `Misc.make_data_dict` and calls `Misc.fragment_data_dict`; every chunk holds the SNPs of one chromosome inside one window of `chunk_size` base pairs. SNP IDs must look like `chromosome_position`.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.fragment_snp_data(path, chunk_size, *, include_snp_ids=False)`.
- **Input:** `path`, `chunk_size` in base pairs, `include_snp_ids`.
- **Output:** `num_snps`, `num_chunks`, and `chunks` (index, number of SNPs, optional IDs).
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import fragment_snp_data

result = fragment_snp_data("snps.txt", 1000, include_snp_ids=True, verbose=False)
print(result["num_snps"], result["num_chunks"])
print([chunk["snp_ids"] for chunk in result["chunks"]])
```

- **Example output:**

```text
3 2
[['g1_chr1_100', 'g1_chr1_200'], ['g2_chr2_50']]
```

## 11) `dnakit.population.expected_sfs_extrapolated()`

- **Purpose:** Expected spectrum extrapolated to an infinitely fine grid.
- **Method:** Evaluates a demographic model on each grid size and extrapolates: `make_extrap_func` for `log_extrapolation=False` and `make_extrap_log_func` (extrapolating the logarithm) for `True`. The model family (one, two or three populations) follows the number of sample sizes. One function covers two work-list rows.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.expected_sfs_extrapolated(model, params, sample_sizes, *, grid_points=(40, 50, 60), log_extrapolation=False)`.
- **Input:** `model`, `params`, `sample_sizes`, `grid_points` (three or more), `log_extrapolation`.
- **Output:** `spectrum`, `model`, `grid_points`, `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import expected_sfs_extrapolated

result = expected_sfs_extrapolated("two_epoch", [2.0, 0.1], [6], grid_points=[30, 40, 50], verbose=False)
print([None if v is None else round(v, 4) for v in result["spectrum"]])
```

- **Example output:**

```text
[None, 1.1253, 0.506, 0.3331, 0.2497, 0.1999, None]
```

## 12) `dnakit.population.log_likelihood_derivatives()`

- **Purpose:** Hessian or gradient of a model's log-likelihood with respect to its parameters.
- **Method:** Calls `Godambe.get_hess` or `Godambe.get_grad` on the log-likelihood of the data under the extrapolated demographic model (multinomial by default, Poisson with `multinom=False`), by finite differences with fractional step `eps`. One function covers two work-list rows.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.log_likelihood_derivatives(data, model, params, *, derivative='hessian', grid_points=(20, 30, 40), eps=0.01, multinom=True)`.
- **Input:** `data`, `model`, `params`, `derivative` (`hessian` or `gradient`), `grid_points`, `eps`, `multinom`.
- **Output:** `value`: a matrix for the Hessian, a list for the gradient.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import (
    expected_sfs_extrapolated,
    log_likelihood_derivatives,
)

truth = expected_sfs_extrapolated("two_epoch", [2.0, 0.3], [10], grid_points=[20, 30, 40], verbose=False)
data = [0.0 if v is None else 500 * v for v in truth["spectrum"]]
grad = log_likelihood_derivatives(data, "two_epoch", [2.0, 0.3], derivative="gradient", verbose=False)
hess = log_likelihood_derivatives(data, "two_epoch", [2.0, 0.3], verbose=False)
print([round(v, 2) for v in grad["value"]])
print([[round(v, 1) for v in row] for row in hess["value"]])
```

- **Example output:**

```text
[0.0, 0.0]
[[-16.7, -30.1], [-30.1, -98.3]]
```

## 13) `dnakit.population.integrate_diffusion_phi()`

- **Purpose:** Integrate the diffusion density phi over time for one to five populations.
- **Method:** Starts from the equilibrium density of an ancestral population (`PhiManip.phi_1D`, or `phi_1D_X` for an X-linked locus), splits it at time zero into the requested number of populations, and integrates it with `Integration.one_pop`, `two_pops`, `three_pops`, `four_pops`, `five_pops` or `one_pop_X` for `duration` (units of 2N generations) with constant sizes, selection, dominance and a migration matrix. With `sample_sizes` the expected spectrum is computed (`Spectrum.from_phi`). One function covers six work-list rows.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.integrate_diffusion_phi(populations, duration, *, nu=None, migration=None, gamma=None, h=None, theta0=1.0, grid_points=30, x_linked=False, breeding_ratio=1.0, mutation_ratio=1.0, sample_sizes=None, return_phi=False)`.
- **Input:** `populations` (1 to 5), `duration`, optional `nu`, `migration` (matrix, rate into i from j), `gamma`, `h`, `theta0`, `grid_points`, `x_linked`, `breeding_ratio`, `mutation_ratio`, `sample_sizes`, `return_phi`.
- **Output:** `phi_shape`, `phi_total`, `grid`, optional `phi` and `spectrum`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import integrate_diffusion_phi

one = integrate_diffusion_phi(1, 0.1, nu=[2.0], grid_points=20, sample_sizes=[6], verbose=False)
print(one["phi_shape"], [None if v is None else round(v, 3) for v in one["spectrum"]])
two = integrate_diffusion_phi(
    2, 0.1, nu=[2.0, 1.0], migration=[[0, 1.0], [0.5, 0]], grid_points=20, sample_sizes=[4, 4], verbose=False
)
print(two["phi_shape"], round(two["spectrum"][1][1], 3))
```

- **Example output:**

```text
[20] [None, 1.253, 0.545, 0.349, 0.257, 0.203, None]
[20, 20] 0.234
```

## 14) `dnakit.population.simulate_low_coverage_reads()`

- **Purpose:** Simulate read counts across a genotype partition.
- **Method:** Calls `LowPass.simulate_reads`: each individual's depth is drawn from the coverage distribution of its population, and heterozygotes split their reads binomially between the alleles. `seed` makes the draw reproducible.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.simulate_low_coverage_reads(coverage_distribution, genotype_partition, sequenced_haplotypes, number_simulations, *, seed=None)`.
- **Input:** `coverage_distribution` (population to probabilities of depths 0, 1, 2, ...), `genotype_partition` (0, 1 or 2 per individual), `sequenced_haplotypes`, `number_simulations`, `seed`.
- **Output:** `reference_reads` and `alternate_reads`: rows are simulated loci, columns individuals.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import simulate_low_coverage_reads

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = simulate_low_coverage_reads(coverage, [0, 1, 2, 1], [8], 4, seed=2, verbose=False)
print(result["reference_reads"])
print(result["alternate_reads"])
```

- **Example output:**

```text
[[3, 0, 0, 2], [3, 0, 0, 3], [2, 1, 0, 3], [1, 2, 0, 3]]
[[0, 0, 3, 1], [0, 2, 2, 1], [0, 1, 4, 0], [0, 1, 2, 1]]
```

## 15) `dnakit.population.subsample_genotypes()`

- **Purpose:** Subsample genotype calls to fewer haplotypes.
- **Method:** Calls `LowPass.subsample_genotypes_1D`: loci with fewer than `subsample_haplotypes / 2` called individuals are dropped, and that many called individuals are drawn at random from every other locus (99 marks a missing call).
- **API:** `dnakit.population.dadi_p22_spectrum_functions.subsample_genotypes(genotype_calls, subsample_haplotypes, *, seed=None)`.
- **Input:** `genotype_calls` (rows = loci, columns = individuals; 0, 1, 2, 99), `subsample_haplotypes`, `seed`.
- **Output:** `genotypes`, `num_loci`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import subsample_genotypes

calls = [[0, 1, 2, 99], [1, 1, 0, 2], [99, 99, 99, 1]]
result = subsample_genotypes(calls, 4, seed=1, verbose=False)
print(result["num_loci"], result["genotypes"])
```

- **Example output:**

```text
2 [[0, 1], [2, 0]]
```

## 16) `dnakit.population.simulate_gatk_calling_spectrum()`

- **Purpose:** Spectrum produced by GATK multi-sample calling of alleles at given counts.
- **Method:** Calls `LowPass.simulate_GATK_multisample_calling`: loci are simulated with the read depths of the coverage distribution, called as GATK would (a variant needs at least two alternate reads), and subsampled. Entry 0 includes the loci that were not called.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.simulate_gatk_calling_spectrum(coverage_distribution, allele_counts, sequenced_haplotypes, subsample_haplotypes, number_simulations, *, inbreeding=None, seed=None)`.
- **Input:** `coverage_distribution`, `allele_counts`, `sequenced_haplotypes`, `subsample_haplotypes` (one per population), `number_simulations`, `inbreeding`, `seed`.
- **Output:** `probabilities`: probability of each observed derived-allele count.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import simulate_gatk_calling_spectrum

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = simulate_gatk_calling_spectrum(coverage, [3], [6], [4], 200, seed=1, verbose=False)
print(result["probabilities"], round(sum(result["probabilities"]), 6))
```

- **Example output:**

```text
[0.095, 0.26, 0.45, 0.17, 0.025] 1.0
```

## 17) `dnakit.population.gatk_no_call_probability()`

- **Purpose:** Probability that GATK multi-sample calling makes no genotype call.
- **Method:** Calls `LowPass.probability_of_no_call_1D_GATK_multisample` for every derived-haplotype count 0 to n.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.gatk_no_call_probability(coverage_distribution, sequenced_haplotypes, *, inbreeding=0.0)`.
- **Input:** `coverage_distribution` (probabilities of depths 0, 1, ...), `sequenced_haplotypes`, `inbreeding`.
- **Output:** `probabilities`, entry k for k derived haplotypes.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import gatk_no_call_probability

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = gatk_no_call_probability(coverage["A"], 6, verbose=False)
print([round(v, 4) for v in result["probabilities"]])
```

- **Example output:**

```text
[1.0, 0.5148, 0.1618, 0.0453, 0.0122, 0.0033, 0.0009]
```

## 18) `dnakit.population.probability_enough_individuals_covered()`

- **Purpose:** Probability that enough individuals are covered to obtain the subsample.
- **Method:** Calls `LowPass.probability_enough_individuals_covered`: for a site already called as variant, the chance that at least `subsample_haplotypes / 2` of the `sequenced_haplotypes / 2` individuals have reads.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.probability_enough_individuals_covered(coverage_distribution, sequenced_haplotypes, subsample_haplotypes)`.
- **Input:** `coverage_distribution`, `sequenced_haplotypes`, `subsample_haplotypes`.
- **Output:** `probability`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import probability_enough_individuals_covered

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = probability_enough_individuals_covered(coverage["A"], 6, 4, verbose=False)
print(round(result["probability"], 4))
```

- **Example output:**

```text
0.9975
```

## 19) `dnakit.population.expected_sfs_low_pass()`

- **Purpose:** Expected spectrum of a demographic model distorted by low-coverage GATK multi-sample calling.
- **Method:** Wraps the extrapolated model with `LowPass.make_low_pass_func_GATK_multisample` and evaluates it: the spectrum for the sequenced haplotypes is transformed to the subsampled size by the probability of no call, genotype-calling errors and projection, with simulations for the entries below `simulation_threshold`.
- **API:** `dnakit.population.dadi_p22_spectrum_functions.expected_sfs_low_pass(model, params, coverage_distribution, sequenced_haplotypes, subsample_haplotypes, *, grid_points=(30, 40, 50), inbreeding=None, simulation_threshold=0.01, simulations_per_frequency=1000, seed=None)`.
- **Input:** `model`, `params`, `coverage_distribution`, `sequenced_haplotypes`, `subsample_haplotypes`, `grid_points`, `inbreeding`, `simulation_threshold`, `simulations_per_frequency`, `seed`.
- **Output:** `spectrum` for the subsampled size and `sample_sizes`.
- **Example code:**

```python
from dnakit.population.dadi_p22_spectrum_functions import expected_sfs_low_pass

coverage = {"A": [0.05, 0.1, 0.2, 0.25, 0.2, 0.1, 0.05, 0.05]}
result = expected_sfs_low_pass(
    "two_epoch", [2.0, 0.1], coverage, [6], [4], grid_points=[20, 30, 40],
    simulations_per_frequency=20, seed=1, verbose=False,
)
print(result["sample_sizes"], [None if v is None else round(v, 3) for v in result["spectrum"]])
```

- **Example output:**

```text
[4] [1.023, 0.48, 0.474, 0.31, 0.115]
```
