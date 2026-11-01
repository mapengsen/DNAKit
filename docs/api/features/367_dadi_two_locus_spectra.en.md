# dadi two-locus spectra and spectrum files

Functions of the `dnakit.population` area that wrap `dadi.TwoLocus` operations on spectra that are already computed (transformations, marginals, mean r^2, LD per bin, log-likelihoods, scaling, conditioning, density sampling) and the text file format of triallelic and two-locus spectra. A two-locus spectrum is a cube `spectrum[AB][Ab][aB]` of haplotype counts for n chromosomes (ab is n minus the three), as nested lists with `None` for masked entries. Requires `dnakit[bio-dadi]` (dadi 2.4.4 was tested); importing `dadi.TwoLocus` makes dadi create `~/.dadi/TwoLocus_cache`. dadi's two-locus demographic models, forward integration and the three optimizers that call them are not offered (they raise `AttributeError` for `dadi.tridiag`, return NaN once that is patched, or call a removed `fmin_l_bfgs_b` argument), and neither are its genotype-spectrum functions, which were written for Python 2.

## 1) `dnakit.population.two_locus_spectrum_transform()`

- **Purpose:** Mask, fold, unfold, project or misidentify a two-locus spectrum.
- **Method:** Calls the dadi routine named by `operation`: `mask_infeasible`, `unfold` and `fold` (`TLSpectrum` methods), `fold_ancestral`, `fold_lr`, `project` and `misidentification` (`TwoLocus.numerics`), or `array_to_spectrum` (apply dadi's infeasible-entry mask to raw values). The input is not modified. One function covers eight work-list rows.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.two_locus_spectrum_transform(spectrum, operation, *, folded=False, sample_size=None, p_misid=None, apply_infeasible_mask=False)`.
- **Input:** A two-locus spectrum `spectrum[AB][Ab][aB]` of side n + 1 (None for masked), `operation`, and as needed `folded`, `sample_size` (project), `p_misid` (misidentification).
- **Output:** The new `spectrum`, `sample_size`, `folded`, `transform`.
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
folded = two_locus_spectrum_transform(spectrum, "fold", verbose=False)
projected = two_locus_spectrum_transform(spectrum, "project", sample_size=4, verbose=False)
print(folded["folded"], folded["sample_size"], projected["sample_size"])
print(round(projected["spectrum"][1][1][1], 4))
```

- **Example output:**

```text
True 6 4
3.48
```

## 2) `dnakit.population.two_locus_marginal_spectra()`

- **Purpose:** Single-locus spectra of locus A and locus B from a two-locus spectrum.
- **Method:** `method='marginal'` calls `TLSpectrum.marginalA` and `marginalB`; `method='single_locus'` calls `TwoLocus.numerics.to_single_locus`, which skips masked entries. One function covers three work-list rows.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.two_locus_marginal_spectra(spectrum, *, method='marginal', folded=False)`.
- **Input:** `spectrum`, `method`, `folded`.
- **Output:** `locus_a` and `locus_b` (entry k = sites with k derived alleles), `sample_size`.
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_marginal_spectra,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
result = two_locus_marginal_spectra(spectrum, verbose=False)
print([None if v is None else round(v, 2) for v in result["locus_a"]])
```

- **Example output:**

```text
[None, 13.5, 18.2, 20.3, 19.5, 15.5, None]
```

## 3) `dnakit.population.two_locus_mean_r2()`

- **Purpose:** Mean r^2 of a two-locus spectrum.
- **Method:** `method='numerics'` calls `TwoLocus.numerics.mean_r2` (a loop over unmasked entries); `method='spectrum_method'` calls `TLSpectrum.mean_r2`. Both weight r^2 = D^2 / (pA (1 - pA) pB (1 - pB)) by the spectrum.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.two_locus_mean_r2(spectrum, *, method='numerics', folded=False)`.
- **Input:** `spectrum`, `method`, `folded`.
- **Output:** `mean_r2` and `sample_size`.
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_mean_r2,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
print(round(two_locus_mean_r2(spectrum, verbose=False)["mean_r2"], 4))
print(round(two_locus_mean_r2(spectrum, method="spectrum_method", verbose=False)["mean_r2"], 4))
```

- **Example output:**

```text
0.3539
0.3539
```

## 4) `dnakit.population.two_locus_ld_per_bin()`

- **Purpose:** Linkage disequilibrium D and r^2 for every haplotype-count bin.
- **Method:** Calls `TwoLocus.numerics.LD_per_bin(n)`; D = pAB - pA pB and r^2 = D^2 / (pA (1 - pA) pB (1 - pB)) with p the haplotype counts divided by n.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.two_locus_ld_per_bin(sample_size)`.
- **Input:** `sample_size` (2 to 40).
- **Output:** `D` and `r2`, cubes of side n + 1 (None where undefined).
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import two_locus_ld_per_bin

result = two_locus_ld_per_bin(4, verbose=False)
print(result["D"][2][0][0], result["r2"][2][0][0])
```

- **Example output:**

```text
0.25 1.0
```

## 5) `dnakit.population.two_locus_log_likelihood()`

- **Purpose:** Poisson or multinomial log-likelihood of a two-locus data spectrum.
- **Method:** Calls `TwoLocus.inference.ll` (Poisson), `ll_multinom` (model rescaled by its optimal scalar), or with `per_bin=True` `ll_per_bin` and `ll_multinom_per_bin`. One function covers four work-list rows.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.two_locus_log_likelihood(model, data, *, multinom=False, per_bin=False)`.
- **Input:** `model` and `data` two-locus spectra of the same size, `multinom`, `per_bin`.
- **Output:** `log_likelihood` or, with `per_bin`, `log_likelihood_per_bin`.
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_log_likelihood,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
data = [[[v * 2 if v is not None else None for v in row] for row in plane] for plane in spectrum]
print(round(two_locus_log_likelihood(spectrum, data, verbose=False)["log_likelihood"], 3))
print(round(two_locus_log_likelihood(spectrum, data, multinom=True, verbose=False)["log_likelihood"], 3))
```

- **Example output:**

```text
-122.296
-88.688
```

## 6) `dnakit.population.two_locus_log_likelihood_over_rho_bins()`

- **Purpose:** Total two-locus log-likelihood over recombination-rate bins.
- **Method:** Calls `TwoLocus.inference.ll_over_rho_bins`, or with `multinom=True` `ll_over_rho_bins_multinom`, which sums the per-bin log-likelihoods. One function covers two work-list rows.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.two_locus_log_likelihood_over_rho_bins(models, data, *, multinom=False)`.
- **Input:** `models` and `data`, lists with one two-locus spectrum per rho bin, and `multinom`.
- **Output:** `log_likelihood`, `num_bins`.
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    two_locus_log_likelihood_over_rho_bins,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
total = two_locus_log_likelihood_over_rho_bins([spectrum, spectrum], [spectrum, spectrum], verbose=False)
print(total["num_bins"], round(total["log_likelihood"], 3))
```

- **Example output:**

```text
2 -139.176
```

## 7) `dnakit.population.scale_two_locus_spectrum()`

- **Purpose:** Optimal scalar between a two-locus model and data, and the scaled model.
- **Method:** Calls `TwoLocus.inference.optimal_sfs_scaling` (data total over model total on bins unmasked in both) and `optimally_scaled_sfs`. One function covers two work-list rows.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.scale_two_locus_spectrum(model, data)`.
- **Input:** `model` and `data` spectra of the same size.
- **Output:** `scaling` and `scaled_model`.
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    scale_two_locus_spectrum,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
data = [[[v * 3 if v is not None else None for v in row] for row in plane] for plane in spectrum]
print(round(scale_two_locus_spectrum(spectrum, data, verbose=False)["scaling"], 3))
```

- **Example output:**

```text
3.0
```

## 8) `dnakit.population.condition_two_locus_spectrum()`

- **Purpose:** Condition a two-locus spectrum on the derived count at locus A.
- **Method:** Calls `TwoLocus.numerics.condition(F, i)`: entry `[a][b]` is the spectrum value with a AB and b aB haplotypes where AB + Ab = i.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.condition_two_locus_spectrum(spectrum, derived_count)`.
- **Input:** `spectrum` and `derived_count` (1 to n - 1).
- **Output:** `spectrum` of shape (derived_count + 1) x (n - derived_count + 1).
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import (
    condition_two_locus_spectrum,
    two_locus_spectrum_transform,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(7)] for j in range(7)] for i in range(7)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
result = condition_two_locus_spectrum(spectrum, 2, verbose=False)
print(len(result["spectrum"]), len(result["spectrum"][0]))
print([round(v, 2) for v in result["spectrum"][1]])
```

- **Example output:**

```text
3 5
[1.2, 1.3, 1.4, 1.5, 1.6]
```

## 9) `dnakit.population.sample_two_locus_density()`

- **Purpose:** Expected two-locus spectrum sampled from a density on a grid.
- **Method:** Calls `TwoLocus.numerics.sample_cached(phi, n, x)` with the uniform grid `linspace(0, 1, g + 1)`. dadi's own forward integration of the density is not offered: it raises `AttributeError: module 'dadi' has no attribute 'tridiag'` for every input and, with the name supplied, returns NaN.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.sample_two_locus_density(phi, sample_size)`.
- **Input:** `phi`, a cube of side g + 1 (at most 41), and `sample_size`.
- **Output:** `spectrum` (cube of side n + 1), `grid_points`, `sample_size`.
- **Example code:**

```python
from dnakit.population.dadi_p22_twolocus_functions import sample_two_locus_density

phi = [[[1.0] * 6 for _ in range(6)] for _ in range(6)]
result = sample_two_locus_density(phi, 4, verbose=False)
print(result["grid_points"], result["sample_size"], round(result["spectrum"][1][1][1], 4))
```

- **Example output:**

```text
5 4 -2.7972
```

## 10) `dnakit.population.write_multilocus_spectrum_file()`

- **Purpose:** Write a triallelic or two-locus spectrum to a dadi text file.
- **Method:** Calls `TriSpectrum.to_file` or `TLSpectrum.to_file` according to `kind`. The file has comment lines, a header with size, folding and extrapolation values, the data line and the mask line (1 = masked). An existing file is never overwritten. One function covers four work-list rows with `read_multilocus_spectrum_file`.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.write_multilocus_spectrum_file(spectrum, destination, *, kind='triallelic', folded_major=False, folded_ancestral=False, extrap_x=None, extrap_t=None, precision=16, comment_lines=())`.
- **Input:** `spectrum`, a new `destination` path, `kind` (`triallelic` or `two_locus`), folding flags, `extrap_x`, `extrap_t`, `precision` (1 to 17), `comment_lines`.
- **Output:** `destination`, `bytes`, `kind`, `sample_size`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population.dadi_p22_dfe_triallele_functions import triallelic_neutral_equilibrium
from dnakit.population.dadi_p22_twolocus_functions import (
    read_multilocus_spectrum_file,
    write_multilocus_spectrum_file,
)

spectrum = triallelic_neutral_equilibrium(20, 6, verbose=False)["spectrum"]
path = Path(tempfile.mkdtemp()) / "tri.fs"
written = write_multilocus_spectrum_file(spectrum, str(path), comment_lines=["neutral"], verbose=False)
print(written["kind"], written["sample_size"])
print(path.read_text().splitlines()[:2])
back = read_multilocus_spectrum_file(str(path), verbose=False)
print(back["comments"], back["folded_major"], back["spectrum"] == spectrum)
```

- **Example output:**

```text
triallelic 6
['# neutral', '6 unfolded_major unfolded_ancestral None None']
['neutral'] False False
```

## 11) `dnakit.population.read_multilocus_spectrum_file()`

- **Purpose:** Read a triallelic or two-locus spectrum from a dadi text file.
- **Method:** Calls `TriSpectrum.from_file` or `TLSpectrum.from_file` according to `kind` (which must match the file). `mask_infeasible` masks the entries that cannot occur.
- **API:** `dnakit.population.dadi_p22_twolocus_functions.read_multilocus_spectrum_file(path, *, kind='triallelic', mask_infeasible=True)`.
- **Input:** `path`, `kind`, `mask_infeasible`.
- **Output:** `spectrum`, `sample_size`, folding flags, `extrap_x`, `extrap_t`, `comments`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population.dadi_p22_twolocus_functions import (
    read_multilocus_spectrum_file,
    two_locus_spectrum_transform,
    write_multilocus_spectrum_file,
)

cube = [[[1.0 + 0.1 * (i + j + k) for k in range(5)] for j in range(5)] for i in range(5)]
spectrum = two_locus_spectrum_transform(cube, "mask_infeasible", verbose=False)["spectrum"]
path = Path(tempfile.mkdtemp()) / "tl.fs"
write_multilocus_spectrum_file(spectrum, str(path), kind="two_locus", extrap_x=0.01, verbose=False)
back = read_multilocus_spectrum_file(str(path), kind="two_locus", verbose=False)
print(back["sample_size"], back["folded"], back["extrap_x"], back["extrap_t"])
print(round(back["spectrum"][1][1][1], 3))
```

- **Example output:**

```text
4 False 0.01 None
1.3
```
