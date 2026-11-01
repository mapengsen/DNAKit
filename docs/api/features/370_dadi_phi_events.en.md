# dadi population densities: admixture and population selection

Population-density (phi) steps of dadi's diffusion pipeline, as Python functions on JSON-friendly nested lists. phi is stored on dadi's default frequency grid with `grid_points` points per axis. Requires `dnakit[bio-dadi]`.

## 1) `dnakit.population.phi_admix_into_population()`

- **Purpose:** Replace one population of a density by a mixture of the other populations.
- **Method:** Pick the dadi routine from the number of populations and the target: `phi_2D_admix_1_into_2` / `phi_2D_admix_2_into_1`, `phi_3D_admix_*_into_*`, `phi_4D_admix_into_1` to `_4`, `phi_5D_admix_into_1` to `_5`. Proportions are given for the other populations in ascending order. Fourteen work-list rows.
- **API:** `dnakit.population.dadi_p23_phi_functions.phi_admix_into_population(phi, *, grid_points, target_population, proportions, output_path=None)`.
- **Input:** `phi` (2 to 5 populations), `grid_points`, 1-based `target_population`, `proportions` (n-1 numbers), optional `output_path`.
- **Output:** The admixed `phi` (same shape), `target_population`, `proportions`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population, phi_admix_into_population

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
mixed = phi_admix_into_population(two["phi"], grid_points=20, target_population=2, proportions=[0.2])
print(mixed["shape"], mixed["target_population"], round(mixed["phi_sum"], 3))
```

- **Example output:**

```text
[20, 20] 2 3554632.261
```

## 2) `dnakit.population.phi_select_populations()`

- **Purpose:** Remove, keep or reorder the populations of a density.
- **Method:** `remove` calls `PhiManip.remove_pop`, `keep` calls `filter_pops`, `reorder` calls `reorder_pops`. Three work-list rows.
- **API:** `dnakit.population.dadi_p23_phi_functions.phi_select_populations(phi, *, grid_points, operation, populations, output_path=None)`.
- **Input:** `phi`, `grid_points`, `operation`, 1-based `populations`, optional `output_path`.
- **Output:** The new `phi`, `num_populations`, `shape`, `select`, `populations`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, phi_split_population, phi_select_populations

one = phi_one_population(grid_points=20, model="neutral")
two = phi_split_population(one["phi"], grid_points=20, event="1_to_2")
for op, pops in (("remove", [2]), ("keep", [2]), ("reorder", [2, 1])):
    out = phi_select_populations(two["phi"], grid_points=20, operation=op, populations=pops)
    print(op, out["shape"])
```

- **Example output:**

```text
remove [20]
keep [20]
reorder [20, 20]
```
