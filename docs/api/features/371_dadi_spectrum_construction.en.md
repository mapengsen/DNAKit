# dadi spectrum construction

Population-density (phi) steps of dadi's diffusion pipeline, as Python functions on JSON-friendly nested lists. phi is stored on dadi's default frequency grid with `grid_points` points per axis. Requires `dnakit[bio-dadi]`.

## 1) `dnakit.population.spectrum_from_phi()`

- **Purpose:** Sample an allele-frequency spectrum from a population density.
- **Method:** Call `Spectrum.from_phi`; when `inbreeding_coefficients` and `ploidies` are given call `Spectrum.from_phi_inbreeding` instead. Two work-list rows.
- **API:** `dnakit.population.dadi_p23_phi_functions.spectrum_from_phi(phi, *, grid_points, sample_sizes, inbreeding_coefficients=None, ploidies=None, mask_corners=True, pop_ids=None, admix_props=None, het_ascertained=None, force_direct=False)`.
- **Input:** `phi`, `grid_points`, `sample_sizes` (one per population), optional inbreeding coefficients and ploidies, `mask_corners`, `pop_ids`, `admix_props`, `het_ascertained`, `force_direct`.
- **Output:** `spectrum` (nested list, masked entries None), `sample_sizes`, `segregating_sites`, `mode`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import phi_one_population, spectrum_from_phi

phi = phi_one_population(grid_points=30, model="neutral")["phi"]
plain = spectrum_from_phi(phi, grid_points=30, sample_sizes=[6])
inbred = spectrum_from_phi(phi, grid_points=30, sample_sizes=[6], inbreeding_coefficients=[0.2], ploidies=[2])
print(plain["mode"], [round(x, 3) for x in plain["spectrum"][1:-1]])
print(inbred["mode"], [round(x, 3) for x in inbred["spectrum"][1:-1]])
```

- **Example output:**

```text
from_phi [1.039, 0.514, 0.339, 0.253, 0.201]
from_phi_inbreeding [0.881, 0.512, 0.327, 0.249, 0.194]
```

## 2) `dnakit.population.scramble_spectrum_pop_ids()`

- **Purpose:** Randomly reassign sampled individuals among the populations of a spectrum.
- **Method:** Call `Spectrum.scramble_pop_ids` under a fixed random seed (global generators restored afterwards).
- **API:** `dnakit.population.dadi_p23_phi_functions.scramble_spectrum_pop_ids(spectrum, *, mask_corners=True, seed=0)`.
- **Input:** `spectrum` with at least two populations, `mask_corners`, `seed`.
- **Output:** `spectrum`, `sample_sizes`, `segregating_sites`, `seed`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import scramble_spectrum_pop_ids

spectrum = [[None, 5.0, 3.0, 1.0], [4.0, 3.0, 2.0, 1.0], [2.0, 2.0, 1.0, 1.0], [1.0, 1.0, 1.0, None]]
first = scramble_spectrum_pop_ids(spectrum, seed=1)
second = scramble_spectrum_pop_ids(spectrum, seed=1)
print(first["sample_sizes"], round(first["segregating_sites"], 3), first["spectrum"] == second["spectrum"])
```

- **Example output:**

```text
[3, 3] 28.0 True
```

## 3) `dnakit.population.spectrum_from_snp_data()`

**Example file:** `toy.snp`

```text
Context Out Allele1 A B Allele2 A B ID
ACA ACA C 4 3 T 2 3 s1
CCG CCG C 5 2 T 1 4 s2
GCT GCT G 3 3 A 3 3 s3
TCA TCA C 6 1 T 0 5 s4
ACG ACG C 2 2 T 4 4 s5
```

- **Purpose:** Build an allele-frequency spectrum from a dadi SNP file, optionally correcting ancestral misidentification.
- **Method:** Read the file with `Misc.make_data_dict`, then call `Spectrum.from_data_dict`; with `misidentification_path` call `Spectrum.from_data_dict_corrected` (Hernandez et al. 2007). Two work-list rows.
- **API:** `dnakit.population.dadi_p23_phi_functions.spectrum_from_snp_data(path, *, pop_ids, projections, polarized=True, mask_corners=True, misidentification_path=None, force_positive=True)`.
- **Input:** SNP file path, `pop_ids`, `projections`, `polarized`, `mask_corners`, optional misidentification file, `force_positive`.
- **Output:** `spectrum`, `sample_sizes`, `segregating_sites`, `folded`, `mode`.
- **Example code:**

```python
from dnakit.population.dadi_p23_phi_functions import spectrum_from_snp_data

result = spectrum_from_snp_data("toy.snp", pop_ids=["A", "B"], projections=[4, 4])
print(result["mode"], result["sample_sizes"], round(result["segregating_sites"], 3))
folded = spectrum_from_snp_data("toy.snp", pop_ids=["A"], projections=[6], polarized=False)
print(folded["folded"], round(folded["segregating_sites"], 3))
```

- **Example output:**

```text
from_data_dict [4, 4] 3.996
True 4.0
```
