# dadi: demographic inference and spectrum plots

This group of functions in `dnakit.population` is a set of adapters over [dadi](https://github.com/RyanGutenkunst/dadi) (diffusion approximations for demographic inference and the distribution of fitness effects) and, for one 3D plot, [MayaVi](https://github.com/enthought/mayavi). Demes graphs are read through [demes](https://github.com/popsim-consortium/demes-python), and the figures are drawn with [matplotlib](https://matplotlib.org). The functions plot one-, two- and three-population site frequency spectra, triallelic and two-locus spectra and bivariate DFEs to PNG files, compute expected spectra of DFE mixture models, fit the parameters of a demes graph to a spectrum, and compute depth-of-coverage distributions. They do not read VCF files or build spectra from sequence data; the spectra are passed in as nested lists, with `None` for masked entries.

dadi is distributed under the 0BSD licence (from the installed package metadata). The licence of MayaVi was not checked here.

The dependencies are optional and are imported only when the matching function is called. If a package is missing, the function raises an error that includes the install command. The examples import from `dnakit.population` and call each function with `verbose=False`. Plot functions write a new PNG file and never overwrite an existing one; the examples write to a temporary directory and print only the checkable facts (sizes, fields, whether the file starts with the PNG signature).

Four of the dadi plotting and DFE functions fail with the installed matplotlib 3.11.2 and dadi 2.4.4 (`plot_trispectrum`, `plot_trispectrum_comparison`, `plot_two_locus_spectrum`, `fit_demes_graph_parameters`) and the variant `mixture_point_pos` of `dfe_mixture_spectrum` fails as well. Each of them has a working counterpart on this page (`plot_triallelic_spectrum`, `plot_triallelic_spectrum_comparison`, `plot_two_locus_spectrum_3d`, `optimize_demes_graph_parameters` and `dfe_mixture_point_pos_spectrum`), and the examples of the four failing functions show the real error.

## Example data

The two demes functions read a demes graph and a dadi inference-options file. All other examples pass their inputs directly.

**Example file:** `graph.yaml`

```text
description: split with two demes
time_units: generations
demes:
  - name: anc
    epochs:
      - start_size: 1000
        end_time: 500
  - name: A
    ancestors: [anc]
    epochs:
      - start_size: 800
        end_time: 0
  - name: B
    ancestors: [anc]
    epochs:
      - start_size: 600
        end_time: 0
```

**Example file:** `options.yaml`

```text
parameters:
  - name: nuA
    values:
      - demes:
          A:
            epochs:
              0: start_size
    lower_bound: 1
    upper_bound: 100000
  - name: nuB
    values:
      - demes:
          B:
            epochs:
              0: start_size
    lower_bound: 1
    upper_bound: 100000
```

## 1) `dnakit.population.coverage_depth_distribution()`

- **Purpose:** Compute the distribution of sequencing depths for each population.
- **Method:** Build the per-site coverage dictionary that dadi's LowPass module expects (one entry per site, holding the depth array of every population) and call `dadi.LowPass.LowPass.compute_cov_dist`. The depths 0 to the largest depth seen and their normalised frequencies come back for each population.
- **API:** `coverage_depth_distribution(site_depths: 'Mapping[str, Sequence[Sequence[int]]]', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib).
- **Input:** `site_depths` maps each population ID (unique, non-empty string) to a list over sites; each site is a non-empty list of non-negative integer depths, one per individual. All populations must have the same non-zero number of sites (at most 1,000,000).
- **Output:** `status`, `operation`, `sites`, `distribution` (for each population the lists `depth` and `frequency`) and `provenance`.
- **Example code:**

```python
from dnakit.population import coverage_depth_distribution

site_depths = {"A": [[1, 2], [2, 2], [3, 1]], "B": [[1, 1], [1, 1], [2, 0]]}
result = coverage_depth_distribution(site_depths, verbose=False)
print(result["status"], result["sites"])
for population, row in result["distribution"].items():
    print(population, row["depth"], [round(f, 3) for f in row["frequency"]])
```

- **Example output:**

```text
completed 3
A [0, 1, 2, 3] [0.0, 0.333, 0.5, 0.167]
B [0, 1, 2] [0.167, 0.667, 0.167]
```

## 2) `dnakit.population.dfe_mixture_spectrum()`

- **Purpose:** Compute the expected two-population spectrum of a DFE (distribution of fitness effects) mixture model.
- **Method:** Build a dadi `Cache1D` for `split_mig_sel_single_gamma` and a `Cache2D` for `split_mig_sel` (gamma bounds 1e-4 to `gamma_max`, `gamma_points` log-spaced values plus the point-mass gammas of the variant), then call the matching mixture function of `dadi.DFE.Cache2D_mod` with a lognormal and a bivariate lognormal, or `dadi.DFE.Vourlaki2022.Vourlaki_mixture` with a gamma and a bivariate independent gamma. Non-finite values in the spectrum are returned as `None`.
- **API:** `dfe_mixture_spectrum(variant: 'str', demographic_params: 'Sequence[float]', sample_sizes: 'Sequence[int]', dfe_params: 'Sequence[float]', theta: 'float' = 1.0, *, grid_points: 'Sequence[int]' = (20, 25, 30), gamma_points: 'int' = 20, gamma_max: 'float' = 2000.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib).
- **Input:** `variant` is `mixture` (4 `dfe_params`), `mixture_symmetric_point_pos` (6), `mixture_point_pos` (8) or `Vourlaki_mixture` (6). `demographic_params` are nu1, nu2, T, m. `sample_sizes` are two integers from 2 to 40, `grid_points` three integers from 8 to 200, `gamma_points` an integer from 5 to 60, `theta` and `gamma_max` are positive. Every number in `demographic_params` and `dfe_params` must be at least 0, so a negative lognormal mu is rejected here (`dfe_mixture_point_pos_spectrum` accepts it). The variant `mixture_point_pos` raises `ConfigurationError` in this function because dadi passes `None` where a number is needed; use `dfe_mixture_point_pos_spectrum` for it. `mixture_symmetric_point_pos` returned an all-`None` spectrum in one trial with these small settings, so check the result for `None`.
- **Output:** `status`, `operation`, `variant`, `demographic_params`, `sample_sizes`, `dfe_params`, `theta`, `grid_points`, `gamma_points`, `gamma_max`, `spectrum` (nested lists, masked entries `None`) and `provenance`. With very coarse grids the spectrum can be numerically unstable, so use grids close to the defaults.
- **Example code:**

```python
import numpy as np

from dnakit.population import dfe_mixture_spectrum

result = dfe_mixture_spectrum(
    "mixture",
    [2.0, 2.0, 0.2, 0.1],
    [4, 4],
    [0.5, 1.0, 0.5, 0.5],
    grid_points=[20, 25, 30],
    gamma_points=10,
    gamma_max=100.0,
    verbose=False,
)
spectrum = np.array([[np.nan if v is None else v for v in row] for row in result["spectrum"]])
print(result["variant"], result["sample_sizes"], spectrum.shape)
print(int(np.isnan(spectrum).sum()), round(float(np.nansum(spectrum)), 2))
```

- **Example output:**

```text
mixture [4, 4] (5, 5)
2 3.05
```

## 3) `dnakit.population.fit_demes_graph_parameters()`

- **Purpose:** Fit the free parameters of a demes graph to an observed frequency spectrum.
- **Method:** Call `dadi.Demes.Inference.optimize` with the graph file, the inference-options file, the spectrum and the integration grid, using the multinomial likelihood (no `uL`), the chosen optimiser and the log parameterisation if `log` is true.
- **API:** `fit_demes_graph_parameters(graph: 'str | Path', inference_options: 'str | Path', sfs: 'Sequence[Any]', pop_ids: 'Sequence[str]', grid_points: 'Sequence[int]', *, method: 'str' = 'fmin', maxiter: 'int' = 1000, log: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). With ruamel.yaml 0.19, which is installed here, dadi's `optimize` reads its options file with `ruamel.yaml.load`, which that release removed, so the call raises `ConfigurationError`. It needs a ruamel.yaml release that still provides `load`; `optimize_demes_graph_parameters` (below) does the same fit and works with ruamel.yaml 0.19.
- **Input:** `graph` and `inference_options` are existing files (a demes YAML graph and a dadi options YAML with a `parameters` list). `sfs` is a 1D, 2D or 3D nested list of non-negative counts (masked entries `None`); `pop_ids` names its axes and must match demes in the graph; `grid_points` is three integers from 8 to 200; `method` is `fmin`, `powell` or `lbfgsb`; `maxiter` is an integer from 1 to 10000; `log` is a bool.
- **Output:** `status`, `operation`, `graph`, `inference_options`, `pop_ids`, `sample_sizes`, `grid_points`, `method`, `maxiter`, `log`, `parameter_names`, `params` (fitted values), `fopt` (objective value, `None` if not finite) and `provenance`.
- **Example code:**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.population import fit_demes_graph_parameters

sfs = [
    [None, 0.5, 0.2, 0.1, None],
    [0.6, 0.9, 0.4, 0.2, 0.1],
    [0.2, 0.3, 0.5, 0.2, 0.1],
    [0.1, 0.2, 0.1, 0.3, 0.1],
    [None, 0.1, 0.1, 0.1, None],
]
try:
    fit_demes_graph_parameters(
        "graph.yaml", "options.yaml", sfs, ["A", "B"], [8, 10, 12], maxiter=3, verbose=False
    )
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi Demes inference failed
```

The example shows the true behaviour in this environment: the fit fails and the error is reported as `ConfigurationError`.

- **Usage with the dependency installed:**

```python
from dnakit.population import fit_demes_graph_parameters

sfs = [
    [None, 0.5, 0.2, 0.1, None],
    [0.6, 0.9, 0.4, 0.2, 0.1],
    [0.2, 0.3, 0.5, 0.2, 0.1],
    [0.1, 0.2, 0.1, 0.3, 0.1],
    [None, 0.1, 0.1, 0.1, None],
]
result = fit_demes_graph_parameters(
    "graph.yaml", "options.yaml", sfs, ["A", "B"], [8, 10, 12], method="powell", maxiter=200
)
print(result["parameter_names"], result["params"], result["fopt"])
```

## 4) `dnakit.population.plot_dfe_bivariate()`

- **Purpose:** Draw a bivariate distribution of fitness effects for negative selection coefficients as a heat map.
- **Method:** Build a logarithmic grid of negative gamma values (the same on both axes) and call `dadi.DFE.Plotting.plot_biv_dfe` with `dadi.DFE.PDFs.biv_lognormal` or `biv_ind_gamma`; the figure is saved as a PNG.
- **API:** `plot_dfe_bivariate(distribution: 'str', params: 'Sequence[float]', image_path: 'str | Path', *, gamma_min: 'float' = 0.001, gamma_max: 'float' = 100.0, grid_size: 'int' = 50, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `distribution` is `biv_lognormal` (`params` = mu, sigma, rho, mu and sigma on the log scale) or `biv_ind_gamma` (`params` = alpha, beta and a third value that is ignored). `gamma_min` and `gamma_max` are positive with `gamma_min` smaller; `grid_size` is an integer from 5 to 200; `colorbar` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `distribution`, `params` and `gamma_grid` (the negative gamma values used on both axes), and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_dfe_bivariate

out = Path(tempfile.mkdtemp()) / "dfe.png"
result = plot_dfe_bivariate("biv_lognormal", [1.0, 1.0, 0.5], out, grid_size=20, verbose=False)
grid = result["gamma_grid"]
print(result["status"], result["distribution"], result["params"])
print(len(grid), round(grid[0], 6), round(grid[-1], 3))
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed biv_lognormal [1.0, 1.0, 0.5]
20 -0.001 -100.0
True
```

## 5) `dnakit.population.plot_dfe_bivariate_point_pos()`

- **Purpose:** Draw a bivariate DFE that also has a positive point mass as a heat map.
- **Method:** Build the same negative gamma grid and call `dadi.DFE.Plotting.plot_biv_point_pos_dfe`, which treats the last four values of `params` as the positive proportion and positive gamma of each population and the earlier values as the parameters of the bivariate distribution.
- **API:** `plot_dfe_bivariate_point_pos(distribution: 'str', params: 'Sequence[float]', image_path: 'str | Path', *, rho: 'float' = 0.0, gamma_min: 'float' = 0.001, gamma_max: 'float' = 100.0, grid_size: 'int' = 50, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `distribution` is `biv_lognormal` or `biv_ind_gamma`. `params` holds at least 5 numbers (distribution parameters followed by the four point-mass values). `rho` (at least -1) links the negative and positive parts; for `biv_lognormal` it is separate from the rho inside `params`. Grid arguments are as in `plot_dfe_bivariate`. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `distribution`, `params`, `rho` and `gamma_grid`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_dfe_bivariate_point_pos

out = Path(tempfile.mkdtemp()) / "dfe_pos.png"
params = [1.0, 1.0, 0.5, 0.1, 5.0, 0.1, 5.0]
result = plot_dfe_bivariate_point_pos(
    "biv_lognormal", params, out, rho=0.5, grid_size=20, colorbar=False, verbose=False
)
print(result["operation"], result["distribution"], result["rho"], len(result["params"]))
print(len(result["gamma_grid"]), round(result["gamma_grid"][-1], 3))
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
plot_dfe_bivariate_point_pos biv_lognormal 0.5 7
20 -100.0
True
```

## 6) `dnakit.population.plot_sfs_1d()`

- **Purpose:** Draw a one-population site frequency spectrum.
- **Method:** Wrap the list in a `dadi.Spectrum` and call `dadi.Plotting.plot_1d_fs`; the figure is saved as a PNG.
- **API:** `plot_sfs_1d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `spectrum` is a non-empty 1D list of non-negative counts; `None` marks a masked entry. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `sample_sizes` and `spectrum` (as dadi holds it, masked entries `None`), and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_1d

out = Path(tempfile.mkdtemp()) / "sfs1d.png"
result = plot_sfs_1d([None, 10, 5, 3, 2, None], out, verbose=False)
print(result["status"], result["sample_sizes"], result["spectrum"])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed [5] [None, 10.0, 5.0, 3.0, 2.0, None]
True
```

## 7) `dnakit.population.plot_sfs_1d_comparison()`

- **Purpose:** Draw a model spectrum against a data spectrum with residuals for one population.
- **Method:** Call `dadi.Plotting.plot_1d_comp_multinom` or `plot_1d_comp_Poisson` with the chosen residual type.
- **API:** `plot_sfs_1d_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, likelihood: 'str' = 'multinom', residual: 'str' = 'Anscombe', plot_masked: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `model` and `data` are 1D lists of non-negative counts. `likelihood` is `multinom` or `poisson`; `residual` is `Anscombe` or `linear`; `plot_masked` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `likelihood`, `residual`, `model` and `data`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_1d_comparison

model = [None, 9, 5, 3, 2, None]
data = [None, 10, 5, 3, 2, None]
folder = Path(tempfile.mkdtemp())
first = plot_sfs_1d_comparison(model, data, folder / "default.png", verbose=False)
second = plot_sfs_1d_comparison(
    model, data, folder / "poisson.png", likelihood="poisson", residual="linear", verbose=False
)
print(first["likelihood"], first["residual"], first["model"])
print(second["likelihood"], second["residual"], second["data"])
print((folder / "default.png").read_bytes()[:4] == b"\x89PNG", (folder / "poisson.png").exists())
```

- **Example output:**

```text
multinom Anscombe [None, 9.0, 5.0, 3.0, 2.0, None]
poisson linear [None, 10.0, 5.0, 3.0, 2.0, None]
True True
```

## 8) `dnakit.population.plot_sfs_2d()`

- **Purpose:** Draw a two-population joint site frequency spectrum as a heat map.
- **Method:** Wrap the list in a `dadi.Spectrum` and call `dadi.Plotting.plot_single_2d_sfs`.
- **API:** `plot_sfs_2d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `spectrum` is a 2D list of non-negative counts (masked entries `None`). `vmin` and `vmax` are optional colour limits with `vmin` smaller; `pop_ids` is two distinct names; `colorbar` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `sample_sizes` and `spectrum`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d

out = Path(tempfile.mkdtemp()) / "sfs2d.png"
spectrum = [[None, 3, 1, 0.5], [2, 1, 0.5, 0.2], [1, 0.5, 0.2, 0.1], [0.5, 0.2, 0.1, None]]
result = plot_sfs_2d(spectrum, out, vmin=0.1, vmax=3.0, pop_ids=["A", "B"], verbose=False)
print(result["status"], result["sample_sizes"])
print(result["spectrum"][0])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed [3, 3]
[None, 3.0, 1.0, 0.5]
True
```

## 9) `dnakit.population.plot_sfs_2d_comparison()`

- **Purpose:** Draw the model, the data and the residual of a two-population spectrum side by side.
- **Method:** Call `dadi.Plotting.plot_2d_comp_multinom` or `plot_2d_comp_Poisson`.
- **API:** `plot_sfs_2d_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, likelihood: 'str' = 'multinom', vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, residual: 'str' = 'Anscombe', adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `model` and `data` are 2D lists of non-negative counts. `likelihood` is `multinom` or `poisson`; `residual` is `Anscombe` or `linear`; `vmin`/`vmax` (with `vmin` smaller), `resid_range` (positive), `pop_ids` and `adjust` (bool) are optional. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `likelihood`, `residual`, `model` and `data`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d_comparison

out = Path(tempfile.mkdtemp()) / "comparison2d.png"
model = [[None, 3, 1, 0.5], [2, 1, 0.5, 0.2], [1, 0.5, 0.2, 0.1], [0.5, 0.2, 0.1, None]]
data = [[None, 2.5, 1, 0.5], [2, 1.2, 0.4, 0.2], [1, 0.6, 0.2, 0.1], [0.5, 0.2, 0.1, None]]
result = plot_sfs_2d_comparison(model, data, out, pop_ids=["A", "B"], verbose=False)
print(result["likelihood"], result["residual"])
print(result["model"][1], result["data"][1])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
multinom Anscombe
[2.0, 1.0, 0.5, 0.2] [2.0, 1.2, 0.4, 0.2]
True
```

## 10) `dnakit.population.plot_sfs_2d_meta_residual()`

- **Purpose:** Draw the meta-residual comparison of a synonymous and a nonsynonymous two-population residual spectrum.
- **Method:** Wrap both residuals in `dadi.Spectrum` objects (negative values allowed) and call `dadi.Plotting.plot_2d_meta_resid`.
- **API:** `plot_sfs_2d_meta_residual(synonymous_residual: 'Sequence[Any]', nonsynonymous_residual: 'Sequence[Any]', image_path: 'str | Path', *, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `synonymous_residual` and `nonsynonymous_residual` are 2D lists of finite numbers (masked entries `None`). `resid_range` is positive if given; `pop_ids` is two distinct names; `adjust` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `synonymous_residual`, `nonsynonymous_residual` and `resid_range`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d_meta_residual

out = Path(tempfile.mkdtemp()) / "meta.png"
synonymous = [[None, 0.5, -0.2, 0.1], [0.1, -0.3, 0.2, 0.0], [-0.1, 0.2, 0.0, -0.1], [0.2, 0.0, -0.1, None]]
nonsynonymous = [[None if v is None else -v for v in row] for row in synonymous]
result = plot_sfs_2d_meta_residual(
    synonymous, nonsynonymous, out, resid_range=1.0, pop_ids=["A", "B"], verbose=False
)
print(result["status"], result["resid_range"])
print(result["synonymous_residual"][0], result["nonsynonymous_residual"][0])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed 1.0
[None, 0.5, -0.2, 0.1] [None, -0.5, 0.2, -0.1]
True
```

## 11) `dnakit.population.plot_sfs_2d_residual()`

- **Purpose:** Draw a two-population residual spectrum (signed values) as a heat map.
- **Method:** Wrap the residual in a `dadi.Spectrum` (negative values allowed) and call `dadi.Plotting.plot_2d_resid`.
- **API:** `plot_sfs_2d_residual(residual: 'Sequence[Any]', image_path: 'str | Path', *, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, colorbar: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `residual` is a 2D list of finite numbers (masked entries `None`). `resid_range` (positive) saturates the colour scale at plus and minus that value; `pop_ids` is two distinct names; `colorbar` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `residual` and `resid_range`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_2d_residual

out = Path(tempfile.mkdtemp()) / "residual2d.png"
residual = [[None, 0.5, -0.2, 0.1], [0.1, -0.3, 0.2, 0.0], [-0.1, 0.2, 0.0, -0.1], [0.2, 0.0, -0.1, None]]
result = plot_sfs_2d_residual(residual, out, resid_range=1.0, pop_ids=["A", "B"], verbose=False)
print(result["status"], result["resid_range"])
print(result["residual"][1])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed 1.0
[0.1, -0.3, 0.2, 0.0]
True
```

## 12) `dnakit.population.plot_sfs_3d()`

- **Purpose:** Draw a three-population spectrum as a logarithmic heat map of slices.
- **Method:** Wrap the list in a `dadi.Spectrum` and call `dadi.Plotting.plot_3d_spectrum`.
- **API:** `plot_sfs_3d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `spectrum` is a 3D list of non-negative counts (masked entries `None`). `vmin` and `vmax` are optional (`vmin` smaller); `pop_ids` is three distinct names. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `sample_sizes` and `spectrum`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_3d

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "sfs3d.png"
result = plot_sfs_3d(cube(6.0), out, vmin=0.1, vmax=6.0, pop_ids=["A", "B", "C"], verbose=False)
print(result["status"], result["sample_sizes"])
print(result["spectrum"][0][0])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed [2, 2, 2]
[None, 3.0, 2.0]
True
```

## 13) `dnakit.population.plot_sfs_3d_comparison()`

- **Purpose:** Draw the model, the data and the residual of a three-population spectrum as slice panels.
- **Method:** Call `dadi.Plotting.plot_3d_comp_multinom` or `plot_3d_comp_Poisson`.
- **API:** `plot_sfs_3d_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, likelihood: 'str' = 'multinom', vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, residual: 'str' = 'Anscombe', adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `model` and `data` are 3D lists of non-negative counts. The other arguments are as in `plot_sfs_2d_comparison`. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `likelihood`, `residual`, `model` and `data`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_3d_comparison

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "comparison3d.png"
result = plot_sfs_3d_comparison(cube(6.0), cube(6.6), out, pop_ids=["A", "B", "C"], verbose=False)
print(result["likelihood"], result["residual"])
print(result["model"][0][0], [None if v is None else round(v, 2) for v in result["data"][0][0]])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
multinom Anscombe
[None, 3.0, 2.0] [None, 3.3, 2.2]
True
```

## 14) `dnakit.population.plot_sfs_3d_pairwise()`

- **Purpose:** Draw the three pairwise two-population marginals of a three-population spectrum.
- **Method:** Call `dadi.Plotting.plot_3d_pairwise`, which marginalises the spectrum over the third population for each pair.
- **API:** `plot_sfs_3d_pairwise(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, adjust: 'bool' = True, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `spectrum` is a 3D list of non-negative counts (masked entries `None`). `vmin`/`vmax` (with `vmin` smaller), `pop_ids` (three names) and `adjust` (bool) are optional. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `sample_sizes` and `spectrum`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_sfs_3d_pairwise

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "pairwise3d.png"
result = plot_sfs_3d_pairwise(cube(6.0), out, pop_ids=["A", "B", "C"], verbose=False)
print(result["status"], result["sample_sizes"])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed [2, 2, 2]
True
```

## 15) `dnakit.population.plot_trispectrum()`

- **Purpose:** Draw a triallelic (three-allele) spectrum as a heat map using dadi's Triallele plotting.
- **Method:** Wrap the list in a `dadi.Triallele` `TriSpectrum` and call `dadi.Triallele.plotting.plot_single_trispectrum`. That upstream function passes a `LogNorm` together with `vmin` and `vmax` to `pcolor`, which matplotlib 3.11 (installed here) rejects, so every call raises `ConfigurationError` in this environment. Use `plot_triallelic_spectrum` (page section 23), which draws the same plot without that call.
- **API:** `plot_trispectrum(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed. It is known to fail with matplotlib 3.11.2, the version in the test environment.
- **Input:** `spectrum` is a 2D list of non-negative counts (masked entries `None`); `folded`, `colorbar` are bools; `vmin`/`vmax` are optional. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** On success `status`, `operation`, `image_path`, `image_format`, `folded`, `spectrum` and `provenance`. Here the call fails with `ConfigurationError`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.population import plot_trispectrum

n = 7
spectrum = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
out = Path(tempfile.mkdtemp()) / "tri.png"
try:
    plot_trispectrum(spectrum, out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
print(out.exists())
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi could not draw plot_trispectrum
False
```

The example shows the true behaviour in this environment: the error is raised and no file is written.

## 16) `dnakit.population.plot_trispectrum_comparison()`

- **Purpose:** Draw two triallelic spectra and their residual side by side using dadi's Triallele plotting.
- **Method:** Wrap both lists in `TriSpectrum` objects and call `dadi.Triallele.plotting.plot_trispectrum_comp`. It draws its two spectrum panels with the same `plot_single_trispectrum` call, so it fails with matplotlib 3.11 in the same way as `plot_trispectrum`. Use `plot_triallelic_spectrum_comparison` (section 24) instead.
- **API:** `plot_trispectrum_comparison(first: 'Sequence[Any]', second: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, title1: 'str' = 'sfs1', title2: 'str' = 'sfs2', title3: 'str' = 'Residual', colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed. It is known to fail with matplotlib 3.11.2, the version in the test environment.
- **Input:** `first` and `second` are 2D lists of non-negative counts; `folded`, `colorbar` are bools; `vmin`, `vmax`, `resid_range` are optional; `title1`, `title2`, `title3` are strings. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** On success `status`, `operation`, `image_path`, `image_format`, `folded`, `first`, `second` and `provenance`. Here the call fails with `ConfigurationError`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.population import plot_trispectrum_comparison

n = 7
spectrum = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
out = Path(tempfile.mkdtemp()) / "tri_comp.png"
try:
    plot_trispectrum_comparison(spectrum, spectrum, out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
print(out.exists())
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi could not draw plot_trispectrum_comparison
False
```

The example shows the true behaviour in this environment: the error is raised and no file is written.

## 17) `dnakit.population.plot_two_locus_comparison()`

- **Purpose:** Draw a two-locus model spectrum against data as 3D panels with residuals.
- **Method:** Wrap both lists in dadi `TLSpectrum` objects and call `dadi.TwoLocus.plotting.plot_3d_comp`.
- **API:** `plot_two_locus_comparison(model: 'Sequence[Any]', data: 'Sequence[Any]', image_path: 'str | Path', *, resid_range: 'float' = 1.0, max_sum: 'int | None' = None, max_index: 'int | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `model` and `data` are 3D lists of non-negative counts. `resid_range` is positive; `max_sum` and `max_index` are optional integers that slice the spectra. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `resid_range`, `model` and `data` (as held by the `TLSpectrum`, which keeps dadi's own mask, so fewer entries than were passed are present), and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_two_locus_comparison

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

def unmasked(cube3):
    return sum(v is not None for plane in cube3 for row in plane for v in row)


out = Path(tempfile.mkdtemp()) / "two_locus_comp.png"
result = plot_two_locus_comparison(cube(6.0), cube(6.6), out, resid_range=1.0, verbose=False)
print(result["status"], result["resid_range"])
print(unmasked(result["model"]), unmasked(result["data"]))
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
completed 1.0
2 2
True
```

## 18) `dnakit.population.plot_two_locus_spectrum()`

- **Purpose:** Draw a two-locus (three-dimensional) spectrum as a 3D scatter plot using dadi's TwoLocus plotting.
- **Method:** Wrap the list in a `TLSpectrum` and call `dadi.TwoLocus.plotting.plot_3d`. That upstream function always ends with `ax.colorbar()`, which a 3D axes does not have, so every call raises `ConfigurationError`. Use `plot_two_locus_spectrum_3d` (section 25), which draws the voxels and adds the colour bar correctly.
- **API:** `plot_two_locus_spectrum(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, max_sum: 'int | None' = None, max_index: 'int | None' = None, colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed. It fails for every input with dadi 2.4.4.
- **Input:** `spectrum` is a 3D list of non-negative counts; `vmin`, `vmax`, `max_sum`, `max_index` are optional; `colorbar` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** On success `status`, `operation`, `image_path`, `image_format`, `spectrum` and `provenance`. Here the call fails with `ConfigurationError`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import ConfigurationError
from dnakit.population import plot_two_locus_spectrum

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "two_locus.png"
try:
    plot_two_locus_spectrum(cube(6.0), out, verbose=False)
except ConfigurationError as error:
    print(type(error).__name__, str(error).split(":")[0])
print(out.exists())
```

- **Example output:**

```text
ConfigurationError [CONFIGURATION_ERROR] dadi could not draw plot_two_locus_spectrum
False
```

The example shows the true behaviour in this environment: the error is raised and no file is written.

## 19) `dnakit.population.plot_3d_spectrum_mayavi()`

- **Purpose:** Draw a three-population spectrum as a MayaVi 3D bar chart saved to a PNG.
- **Method:** Check the inputs, then draw one bar per entry that is not masked and is at least `vmin`, coloured by the logarithm of the count scaled between `vmin` and `vmax`, with `mayavi.mlab.barchart`, and save the figure with `mlab.savefig`. It follows dadi's `plot_3d_spectrum_mayavi`.
- **API:** `plot_3d_spectrum_mayavi(spectrum: 'Sequence[Any]', output_path: 'str', *, vmin: 'float | None' = None, vmax: 'float | None' = None, pop_ids: 'Sequence[str] | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package mayavi (extra `dnakit[bio-mayavi]`, which also installs dadi) and a display or an offscreen VTK build. Without mayavi it raises `BackendUnavailableError` with an install hint; MayaVi is not installed in the test environment.
- **Input:** `spectrum` is a 3D list of non-negative counts with at least one positive unmasked entry. `vmin` and `vmax` are positive with `vmin` smaller and default to the smallest and largest positive count. `pop_ids` is three distinct names (default `pop0`, `pop1`, `pop2`). `output_path` must end with `.png` and must not exist; its parent directory is created.
- **Output:** A dictionary with `path` and `data`; `data` holds `pop_ids`, `sample_sizes`, `vmin`, `vmax`, `plotted_entries`, `spectrum` and `provenance` (this function does not return `status` or `operation`).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.exceptions import BackendUnavailableError
from dnakit.population import plot_3d_spectrum_mayavi

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

out = Path(tempfile.mkdtemp()) / "bars.png"
try:
    plot_3d_spectrum_mayavi(cube(6.0), out, verbose=False)
except BackendUnavailableError as error:
    print(type(error).__name__, str(error).split(":")[0])
```

- **Example output:**

```text
BackendUnavailableError [BACKEND_UNAVAILABLE] mayavi.mlab could not be imported
```

The example shows the true behaviour without MayaVi. If MayaVi is installed, the call draws the figure and nothing is printed.

- **Usage with the dependency installed:**

```python
from dnakit.population import plot_3d_spectrum_mayavi

corners = {(0, 0, 0), (2, 2, 2)}
spectrum = [
    [[None if (i, j, k) in corners else 6.0 / (1 + i + j + k) for k in range(3)] for j in range(3)]
    for i in range(3)
]
result = plot_3d_spectrum_mayavi(spectrum, "bars.png", vmin=0.5, pop_ids=["A", "B", "C"])
print(result["path"], result["data"]["plotted_entries"], result["data"]["vmin"], result["data"]["vmax"])
```

## 20) `dnakit.population.vourlaki_visual_validation_png()`

- **Purpose:** Redraw dadi's Vourlaki (2022) DFE validation figure to a PNG and return the spectra it shows.
- **Method:** Build the `Cache1D` and `Cache2D` caches for the dadi visual-validation model (nu1 = nu2 = 2, T = 0.2, m = 0, gamma bounds 1e-4 to 200, a point mass at gamma 10), compute `Vourlaki_mixture` for four scenarios (panels A to D: ppos_wild, pchange and pchange_pos = 0/0/0, 1/0/0, 0/1/0, 0/1/0.9; alpha = 0.2, beta = 10), plot each two-population spectrum with `plot_single_2d_sfs` and the wild and domesticated marginals on log-scale panels, and save one 10 by 8 inch figure at 150 dpi.
- **API:** `vourlaki_visual_validation_png(output_path: 'str', *, sample_sizes: 'Sequence[int]' = (10, 10), grid_points: 'Sequence[int]' = (30, 35, 40), gamma_points: 'int' = 30, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `output_path` must end with `.png` and must not exist; its parent directory is created. `sample_sizes` are two integers from 2 to 40, `grid_points` three integers from 8 to 200, `gamma_points` an integer from 5 to 60. The defaults match upstream (10, 10), (30, 35, 40) and 30; those take longer than the small values used below.
- **Output:** A dictionary with `path` and `data`; `data` holds `sample_sizes`, `grid_points`, `gamma_points`, `demographic_params`, `alpha`, `beta`, `gamma_pos`, `panels` (for each panel `panel`, `ppos_wild`, `pchange`, `pchange_pos`, `spectrum`, `wild_marginal`, `domesticate_marginal`) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

import numpy as np

from dnakit.population import vourlaki_visual_validation_png

out = Path(tempfile.mkdtemp()) / "figures" / "validation.png"
result = vourlaki_visual_validation_png(
    out, sample_sizes=[4, 4], grid_points=[20, 25, 30], gamma_points=10, verbose=False
)
data = result["data"]
print(sorted(result), out.read_bytes()[:4] == b"\x89PNG")
print(data["sample_sizes"], data["demographic_params"], data["alpha"], data["beta"], data["gamma_pos"])
for panel in data["panels"]:
    spectrum = np.array([[np.nan if v is None else v for v in row] for row in panel["spectrum"]])
    print(panel["panel"], panel["ppos_wild"], panel["pchange"], panel["pchange_pos"],
          spectrum.shape, round(float(np.nansum(spectrum)), 2), len(panel["wild_marginal"]))
```

- **Example output:**

```text
['data', 'path'] True
[4, 4] [2, 2, 0.2, 0] 0.2 10.0 10.0
A 0.0 0.0 0.0 (5, 5) 3.43 5
B 1.0 0.0 0.0 (5, 5) 6.41 5
C 0.0 1.0 0.0 (5, 5) 4.93 5
D 0.0 1.0 0.9 (5, 5) 5.55 5
```

## 21) `dnakit.population.dfe_mixture_point_pos_spectrum()`

- **Purpose:** Compute the expected two-population spectrum of dadi's point-mass DFE mixture (`mixture_point_pos`).
- **Method:** Build `Cache1D` and `Cache2D` caches for `split_mig_sel_single_gamma` and `split_mig_sel` with the two positive gammas as extra cache values, compute the one-dimensional part with `integrate_point_pos` and the two-dimensional part with `integrate_point_pos` of the 2D cache, and combine them as (1 - weight) times the first plus weight times the second. It passes `rho` where upstream `mixture_point_pos` passes `None` and fails.
- **API:** `dfe_mixture_point_pos_spectrum(demographic_params: 'Sequence[float]', sample_sizes: 'Sequence[int]', dfe_params: 'Sequence[float]', theta: 'float' = 1.0, *, grid_points: 'Sequence[int]' = (20, 25, 30), gamma_points: 'int' = 20, gamma_max: 'float' = 2000.0, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib).
- **Input:** `demographic_params` are nu1, nu2, T, m (each at least 0). `dfe_params` are exactly 8 finite numbers: lognormal mu and sigma, rho, positive proportion and gamma for population 1, positive proportion and gamma for population 2, and the weight of the 2D part (between 0 and 1); mu may be negative. The two positive gammas must differ, otherwise the caches have a different shape and dadi raises an error. Other arguments are as in `dfe_mixture_spectrum`.
- **Output:** `status`, `operation`, `variant` (`mixture_point_pos`), `demographic_params`, `sample_sizes`, `dfe_params`, `theta`, `grid_points`, `gamma_points`, `gamma_max`, `spectrum` and `provenance`.
- **Example code:**

```python
import numpy as np

from dnakit.population import dfe_mixture_point_pos_spectrum

dfe_params = [-1.0, 0.6, 0.4, 0.2, 5.0, 0.3, 3.0, 0.6]
result = dfe_mixture_point_pos_spectrum(
    [2.0, 2.0, 0.2, 0.1],
    [4, 4],
    dfe_params,
    grid_points=[20, 25, 30],
    gamma_points=10,
    gamma_max=100.0,
    verbose=False,
)
spectrum = np.array([[np.nan if v is None else v for v in row] for row in result["spectrum"]])
print(result["variant"], result["sample_sizes"], spectrum.shape)
print(int(np.isnan(spectrum).sum()), round(float(np.nansum(spectrum)), 2))
```

- **Example output:**

```text
mixture_point_pos [4, 4] (5, 5)
2 5.22
```

## 22) `dnakit.population.optimize_demes_graph_parameters()`

- **Purpose:** Fit the free parameters of a demes graph to a two- or three-population spectrum.
- **Method:** Run the steps of dadi's `Demes.Inference.optimize` itself (parameter set-up, bounds, constraints, optional perturbation, log parameterisation, `scipy.optimize.fmin`, `fmin_powell` or `fmin_l_bfgs_b`) but read the options file with `ruamel.yaml.YAML(typ="safe")`, so it works with ruamel.yaml 0.19. If `output_path` is given, the fitted graph is written with `demes.dump`.
- **API:** `optimize_demes_graph_parameters(graph_path: 'str | Path', inference_options_path: 'str | Path', sfs: 'Sequence[Any]', pop_ids: 'Sequence[str]', grid_points: 'Sequence[int]', *, maxiter: 'int' = 1000, perturb: 'float' = 0, log: 'bool' = True, method: 'str' = 'fmin', uL: 'float | None' = None, output_path: 'str | Path | None' = None, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib).
- **Input:** `graph_path` and `inference_options_path` are existing files. `sfs` is a 1D to 3D nested list of non-negative counts and `pop_ids` names its axes. `grid_points` is three integers from 8 to 200; `method` is `fmin`, `powell` or `lbfgsb`; `maxiter` is an integer from 1 to 10000; `perturb` is a number of at least 0; `log` is a bool; `uL` is an optional positive number; `output_path` is an optional new file.
- **Output:** `status`, `operation`, `method`, `log`, `maxiter`, `parameter_names`, `parameters` (fitted values), `objective`, `output_path` (path of the fitted YAML or `None`) and `provenance`.
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import optimize_demes_graph_parameters

sfs = [
    [None, 0.5, 0.2, 0.1, None],
    [0.6, 0.9, 0.4, 0.2, 0.1],
    [0.2, 0.3, 0.5, 0.2, 0.1],
    [0.1, 0.2, 0.1, 0.3, 0.1],
    [None, 0.1, 0.1, 0.1, None],
]
fitted = Path(tempfile.mkdtemp()) / "fitted.yaml"
result = optimize_demes_graph_parameters(
    "graph.yaml", "options.yaml", sfs, ["A", "B"], [8, 10, 12],
    maxiter=3, output_path=fitted, verbose=False,
)
print(result["method"], result["log"], result["parameter_names"])
print([round(v) for v in result["parameters"]], round(result["objective"], 2))
text = fitted.read_text()
print(text.splitlines()[0], "start_size: 868" in text)
```

- **Example output:**

```text
fmin True ['nuA', 'nuB']
[800, 869] 11.73
description: split with two demes True
```

## 23) `dnakit.population.plot_triallelic_spectrum()`

- **Purpose:** Draw a triallelic spectrum as a log-scale heat map in a new PNG.
- **Method:** Fold the spectrum first when `folded` is true (as dadi's Triallele `_fold` does), then draw it with matplotlib `pcolor` and a `LogNorm` alone. Entries below `vmin` are not drawn. It replaces dadi's `plot_single_trispectrum`, which fails with matplotlib 3.11.
- **API:** `plot_triallelic_spectrum(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, colorbar: 'bool' = False, cmap: 'str' = 'cubehelix_r', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `spectrum` is a square 2D list of non-negative counts (masked entries `None`); a folded spectrum needs at least 3 by 3, and `vmax` must be greater than `vmin`. `vmin` defaults to the smallest positive count and `vmax` to the largest value; `cmap` is a matplotlib colormap name (default `cubehelix_r`); `colorbar` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `folded`, `vmin`, `vmax`, `cmap`, `spectrum` (the entries that were drawn, the others `None`), and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_triallelic_spectrum

n = 7
first = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]

def drawn(grid):
    return sum(v is not None for row in grid for v in row)


folder = Path(tempfile.mkdtemp())
plain = plot_triallelic_spectrum(first, folder / "tri.png", verbose=False)
folded = plot_triallelic_spectrum(first, folder / "tri_folded.png", folded=True, colorbar=True, verbose=False)
print(plain["folded"], plain["vmin"], plain["vmax"], plain["cmap"], drawn(plain["spectrum"]))
print(folded["folded"], folded["vmin"], folded["vmax"], drawn(folded["spectrum"]))
print((folder / "tri.png").read_bytes()[:4] == b"\x89PNG", (folder / "tri_folded.png").exists())
```

- **Example output:**

```text
False 1.0 22.0 cubehelix_r 36
True 14.0 36.0 6
True True
```

## 24) `dnakit.population.plot_triallelic_spectrum_comparison()`

- **Purpose:** Draw two triallelic spectra and their residual in three stacked panels.
- **Method:** Draw both spectra on one log colour scale, compute the residual (first - second) / sqrt(first), and draw it with the diverging `RdBu` colour map centred on zero. It replaces dadi's `plot_trispectrum_comp`, which fails with matplotlib 3.11. With `folded=True` the call currently raises `ConfigurationError` ("the spectrum is already folded"), because the residual of the already folded spectra is folded a second time.
- **API:** `plot_triallelic_spectrum_comparison(first: 'Sequence[Any]', second: 'Sequence[Any]', image_path: 'str | Path', *, folded: 'bool' = False, vmin: 'float | None' = None, vmax: 'float | None' = None, resid_range: 'float | None' = None, colorbar: 'bool' = False, titles: 'Sequence[str]' = ('sfs1', 'sfs2', 'Residual'), cmap: 'str' = 'cubehelix_r', verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `first` and `second` are square 2D lists of non-negative counts of the same shape. `resid_range` is positive if given (default: the largest absolute residual); `titles` is three strings; `vmin`, `vmax`, `cmap` and `colorbar` are as in `plot_triallelic_spectrum`. Keep `folded=False`. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `folded`, `vmin`, `vmax`, `resid_range`, `first`, `second`, `residual`, and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_triallelic_spectrum_comparison

n = 7
first = [[0.0 if i == 0 or j == 0 else max(1.0, 30.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
second = [[0.0 if i == 0 or j == 0 else max(1.0, 28.0 - 4 * (i + j)) for j in range(n)] for i in range(n)]
out = Path(tempfile.mkdtemp()) / "tri_comp.png"
result = plot_triallelic_spectrum_comparison(
    first, second, out, titles=["model", "data", "Residual"], verbose=False
)
print(result["vmin"], result["vmax"], round(result["resid_range"], 4))
print([round(v, 4) for v in result["residual"][1][1:4]])
print(out.read_bytes()[:4] == b"\x89PNG")
```

- **Example output:**

```text
1.0 22.0 0.8165
[0.4264, 0.4714, 0.5345]
True
```

## 25) `dnakit.population.plot_two_locus_spectrum_3d()`

- **Purpose:** Draw a three-dimensional two-locus spectrum as coloured voxels in a new PNG.
- **Method:** Hide entries below `vmin`, apply `max_sum` and `max_index`, colour each voxel with the scrunched `cubehelix_r` map on a log scale, draw only the faces that touch no other voxel with `Poly3DCollection`, and add the colour bar with `figure.colorbar`. It replaces dadi's `TwoLocus.plotting.plot_3d`, which fails with a 3D axes.
- **API:** `plot_two_locus_spectrum_3d(spectrum: 'Sequence[Any]', image_path: 'str | Path', *, vmin: 'float | None' = None, vmax: 'float | None' = None, max_sum: 'int | None' = None, max_index: 'int | None' = None, colorbar: 'bool' = False, verbose: bool = True) -> 'dict[str, Any]'`.
- **Requirements:** Needs the optional package dadi (extra `dnakit[bio-dadi]`, which also installs demes; dadi itself depends on NumPy, SciPy and matplotlib). Figures are drawn with the matplotlib Agg backend, so no display is needed.
- **Input:** `spectrum` is a cubic 3D list of non-negative counts (masked entries `None`). `vmin` defaults to the smallest positive count and `vmax` to the largest value; `max_sum` hides voxels whose index sum exceeds it; `max_index` hides voxels whose three indices all exceed it (both optional integers of at least 0); `colorbar` is a bool. `image_path` must end with `.png`, must not exist yet, and its parent directory must exist.
- **Output:** `status`, `operation`, `image_path` (absolute path of the new PNG), `image_format` (`png`), `vmin`, `vmax`, `max_sum`, `max_index`, `drawn_entries` (number of voxels drawn), `spectrum` (entries at or above `vmin`), and `provenance` (the dadi version and the drawing convention).
- **Example code:**

```python
import tempfile
from pathlib import Path

from dnakit.population import plot_two_locus_spectrum_3d

def cube(scale):
    corners = {(0, 0, 0), (2, 2, 2)}
    return [
        [[None if (i, j, k) in corners else scale / (1 + i + j + k) for k in range(3)] for j in range(3)]
        for i in range(3)
    ]

folder = Path(tempfile.mkdtemp())
full = plot_two_locus_spectrum_3d(cube(6.0), folder / "voxels.png", colorbar=True, verbose=False)
cut = plot_two_locus_spectrum_3d(cube(6.0), folder / "voxels_cut.png", max_sum=2, verbose=False)
print(full["vmin"], full["vmax"], full["drawn_entries"])
print(cut["max_sum"], cut["drawn_entries"])
print((folder / "voxels.png").read_bytes()[:4] == b"\x89PNG", (folder / "voxels_cut.png").exists())
```

- **Example output:**

```text
1.0 3.0 25
2 9
True True
```
