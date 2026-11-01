# dadi population densities: construct and split

Population-density (phi) steps of dadi's diffusion pipeline, as Python functions on JSON-friendly nested lists. phi is stored on dadi's default frequency grid with `grid_points` points per axis. Requires `dnakit[bio-dadi]`.

## 1) `dnakit.population.phi_one_population()`

- **Purpose:** Build the equilibrium population density phi of one population on dadi's frequency grid.
- **Method:** Pick the dadi routine named by `model` and call it on `Numerics.default_grid(grid_points)`: `dominance` = `PhiManip.phi_1D` (numerical quadrature), `genic` = `phi_1D_genic` (analytic series), `neutral` = `phi_1D_snm`, `x_linked` = `phi_1D_X`. This one function covers four work-list rows.
- **API:** `dnakit.population.dadi_p23_phi_functions.phi_one_population(*, grid_points, model='dominance', nu=1.0, theta0=1.0, gamma=0.0, h=0.5, theta=None, beta=1.0, alpha=1.0, output_path=None)`.
- **Input:** `grid_points` (8 to 400), `model`, size `nu`, mutation `theta0`, selection `gamma`, dominance `h`, optional `theta`, `beta`, `alpha`, optional new `.npy` `output_path`.
- **Output:** `phi` (list of grid_points floats), `phi_sum`, `shape`, `phi_path`, `provenance`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population

result = phi_one_population(grid_points=20, model="genic", gamma=1.0)
print(result["model"], result["shape"], round(result["phi_sum"], 3))
print([round(x, 3) for x in result["phi"][:3]])
```

- **Example output:**

```text
genic [20] 5675.193
[2257.854, 2257.854, 681.0]
```

## 2) `dnakit.population.phi_split_population()`

- **Purpose:** Split one population of a density into two (or three) populations.
- **Method:** Call `PhiManip.phi_1D_to_2D`, `phi_2D_to_3D_split_2` or `phi_2D_to_3D_split_1` according to `event` (`1_to_2`, `2_to_3_split_2`, `2_to_3_split_1`). The new population starts as a copy of its parent. Three work-list rows.
- **API:** `dnakit.population.dadi_p23_phi_functions.phi_split_population(phi, *, grid_points, event, output_path=None)`.
- **Input:** `phi` as nested list or `.npy` path, `grid_points`, `event`, optional `output_path`.
- **Output:** The new `phi` with one more axis, `num_populations`, `shape`, `event`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
three = phi_split_population(two["phi"], grid_points=20, event="2_to_3_split_1")
print(two["shape"], three["shape"], three["event"])
```

- **Example output:**

```text
[20, 20] [20, 20, 20] 2_to_3_split_1
```

## 3) `dnakit.population.phi_create_admixed_population()`

- **Purpose:** Add a new population formed by admixture of all existing populations.
- **Method:** With 2, 3 or 4 populations call `PhiManip.phi_2D_to_3D_admix`, `phi_3D_to_4D` or `phi_4D_to_5D`. `proportions` gives the shares of populations 1 to n-1; the last population supplies the rest. Three work-list rows.
- **API:** `dnakit.population.dadi_p23_phi_functions.phi_create_admixed_population(phi, *, grid_points, proportions, output_path=None)`.
- **Input:** `phi`, `grid_points`, `proportions` (n-1 numbers in [0, 1] summing to at most 1), optional `output_path`.
- **Output:** The new `phi`, `num_populations`, `shape`, `proportions`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population, phi_create_admixed_population

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
three = phi_create_admixed_population(two["phi"], grid_points=20, proportions=[0.3])
print(three["shape"], three["proportions"], round(three["phi_sum"], 3))
```

- **Example output:**

```text
[20, 20, 20] [0.3] 4434397133.35
```
