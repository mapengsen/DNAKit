# dadi DFE-integrated and triallelic spectra

Functions of the `dnakit.population` area that wrap `dadi.DFE` (spectra integrated over a distribution of fitness effects) and `dadi.Triallele` (the neutral triallelic density, its sampling, and the folding, unfolding and masking of triallelic spectra). Spectra are nested lists; an entry that dadi masks is `None`. Requires `dnakit[bio-dadi]` (dadi 2.4.4 was tested). dadi's other triallelic functions fail upstream for every input and are not offered: the demographic models call `np.float`, the forward integration returns NaN or zeros, and ancestral misidentification, projection and optimal scaling raise `NameError` or `AttributeError`.

## 1) `dnakit.population.dfe_integrated_spectrum()`

- **Purpose:** Expected spectrum integrated over a distribution of fitness effects (DFE).
- **Method:** Builds dadi's selection cache (`Cache1D` or `Cache2D`) for a demographic model with selection and integrates it with the chosen method: `negative_1d` = `Cache1D.integrate`, `point_pos_1d` = `integrate_point_pos`, `bivariate_2d` = `Cache2D.integrate`, `bivariate_point_pos_2d` = `integrate_point_pos`, `bivariate_symmetric_point_pos_2d` = `integrate_symmetric_point_pos`. One function covers five work-list rows; the model and the density are chosen by name.
- **API:** `dnakit.population.dadi_p22_dfe_triallele_functions.dfe_integrated_spectrum(kind, model, demographic_params, sample_sizes, dfe_params, *, distribution=None, theta=1.0, n_positive=1, grid_points=(20, 25, 30), gamma_points=20, gamma_max=2000.0, exterior_int=True)`.
- **Input:** `kind`, `model` (a `dadi.DFE.DemogSelModels` name), the demographic parameters, one or two sample sizes, the DFE parameters (the last values are the positive point masses for the point-mass kinds), and optional `distribution` (`gamma`, `lognormal`, `exponential`, `beta`, `biv_lognormal`, `biv_ind_gamma`), `theta`, `grid_points` (three grid sizes), `gamma_points`, `gamma_max`, `exterior_int`.
- **Output:** `spectrum` (nested list, masked entries are None), `kind`, `model`, `distribution`, `sample_sizes`, the cache settings and `provenance`.
- **Example code:**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import dfe_integrated_spectrum

result = dfe_integrated_spectrum(
    "negative_1d", "two_epoch_sel", [2.0, 0.1], [6], [0.2, 100.0],
    grid_points=[30, 40, 50], gamma_points=8, gamma_max=100.0, verbose=False,
)
print(result["kind"], result["distribution"], result["sample_sizes"])
print([None if v is None else round(v, 3) for v in result["spectrum"]])
```

- **Example output:**

```text
negative_1d gamma [6]
[None, 1.123, 0.42, 0.253, 0.179, 0.137, None]
```

## 2) `dnakit.population.triallelic_neutral_equilibrium()`

- **Purpose:** Neutral equilibrium triallelic density and spectrum.
- **Method:** Computes the exact neutral density of Jenkins et al. (2014) on the grid (`Triallele.integration.equilibrium_neutral_exact`) and samples it to `sample_size` chromosomes (`Triallele.numerics.sample`), optionally folding by major and minor derived allele. This is the `Triallele.demographics.equilibrium` spectrum for neutral alleles; a density given as `phi` is sampled directly. dadi's forward integration (`advance`, used for selection and the epoch models) returns NaN or zeros for every input in the installed environment, and the models call the removed `np.float`, so only the neutral equilibrium is offered. One function covers three work-list rows.
- **API:** `dnakit.population.dadi_p22_dfe_triallele_functions.triallelic_neutral_equilibrium(grid_points, sample_size=None, *, folded=False, phi=None, return_density=False)`.
- **Input:** `grid_points` (10 to 120), optional `sample_size`, `folded`, `phi` (square nested list), `return_density`.
- **Output:** `spectrum` (square nested list), `sample_size`, `folded_major`, `folded_ancestral`, `density_total`, optional `phi`.
- **Example code:**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import triallelic_neutral_equilibrium

result = triallelic_neutral_equilibrium(30, 8, verbose=False)
print(result["sample_size"], result["folded_major"], round(result["density_total"], 1))
print([None if v is None else round(v, 3) for v in result["spectrum"][1][:4]])
```

- **Example output:**

```text
8 False 12675.3
[None, 0.758, 0.434, 0.291]
```

## 3) `dnakit.population.alternate_mutation_triallelic_spectrum()`

- **Purpose:** Exact neutral triallelic spectrum when both mutations enter at one site.
- **Method:** Calls `Triallele.integration.alt_mut_mech_sample_spectrum`, the Jenkins et al. (2014) solution for the alternate mutation mechanism. Population size changes do not alter this neutral spectrum.
- **API:** `dnakit.population.dadi_p22_dfe_triallele_functions.alternate_mutation_triallelic_spectrum(sample_size)`.
- **Input:** `sample_size`, the number of sampled chromosomes (3 to 100).
- **Output:** `spectrum` (square nested list, infeasible entries None) and `sample_size`.
- **Example code:**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import alternate_mutation_triallelic_spectrum

result = alternate_mutation_triallelic_spectrum(6, verbose=False)
print(result["spectrum"][1])
```

- **Example output:**

```text
[None, 0.5, 0.125, 0.05, 0.025, None, None]
```

## 4) `dnakit.population.triallelic_spectrum_transform()`

- **Purpose:** Mask, unfold, fold or log-transform a triallelic spectrum.
- **Method:** Builds a `TriSpectrum` and calls the method named by `operation`: `mask_infeasible`, `unfold`, `fold_major` (by major and minor derived allele), `fold_ancestral` (ancestral state unknown, tracking the two minor counts) or `log`. One function covers five work-list rows.
- **API:** `dnakit.population.dadi_p22_dfe_triallele_functions.triallelic_spectrum_transform(spectrum, operation, *, folded_major=False, folded_ancestral=False)`.
- **Input:** `spectrum` (square nested list, None for masked entries), `operation`, and the current folding flags `folded_major`, `folded_ancestral`.
- **Output:** The transformed `spectrum`, `sample_size`, `folded_major`, `folded_ancestral`, `transform`.
- **Example code:**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import (
    triallelic_neutral_equilibrium,
    triallelic_spectrum_transform,
)

spectrum = triallelic_neutral_equilibrium(20, 6, verbose=False)["spectrum"]
folded = triallelic_spectrum_transform(spectrum, "fold_major", verbose=False)
print(folded["folded_major"], folded["folded_ancestral"])
print([None if v is None else round(v, 3) for v in folded["spectrum"][2]])
back = triallelic_spectrum_transform(folded["spectrum"], "unfold", folded_major=True, verbose=False)
print(back["folded_major"])
```

- **Example output:**

```text
True False
[None, 0.851, 0.247, None, None, None, None]
False
```

## 5) `dnakit.population.triallelic_segregating_sites()`

- **Purpose:** Number of triallelic segregating sites.
- **Method:** Sums the unmasked entries of the spectrum (`TriSpectrum.S`). `apply_infeasible_mask=True` first masks the entries that cannot occur, as dadi does when it builds a spectrum.
- **API:** `dnakit.population.dadi_p22_dfe_triallele_functions.triallelic_segregating_sites(spectrum, *, folded_major=False, folded_ancestral=False, apply_infeasible_mask=False)`.
- **Input:** `spectrum` (square nested list, None for masked), `folded_major`, `folded_ancestral`, `apply_infeasible_mask`.
- **Output:** `segregating_sites` and `sample_size`.
- **Example code:**

```python
from dnakit.population.dadi_p22_dfe_triallele_functions import (
    triallelic_neutral_equilibrium,
    triallelic_segregating_sites,
)

spectrum = triallelic_neutral_equilibrium(20, 6, verbose=False)["spectrum"]
result = triallelic_segregating_sites(spectrum, verbose=False)
print(result["sample_size"], round(result["segregating_sites"], 4))
```

- **Example output:**

```text
6 3.1528
```
