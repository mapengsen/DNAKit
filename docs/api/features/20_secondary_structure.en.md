# Secondary structure properties

Resolve DNA secondary structure and calculate base pairing, stem-loop, pairing probability, accessibility, collective defects, and target structure probability.

Native functionality supports dot-bracket parsing and pairing probability derived metrics; NUPACK analysis is only performed if the user has licensed and independently installed NUPACK 4, projects do not automatically download, install or silently call NUPACK.

The paper sources and internal formulas for various calculations can be found in [FAQ: Calculation basis and references for secondary structure properties ](../../faq.md#secondary-structure-references).

## 1) Dot-bracket structure analysis

- **Function:** Parse the dot-bracket string, restore the base pairing relationship and count structural elements such as stem, hairpin, unpaired bases, etc., for checking and summarizing existing secondary structure predictions.
- **Calculation method:** Use stack to pair open and close positions according to bracket type; continuous `(i,j)、(i+1,j−1)` pairings are merged into stems, and the innermost chain interval without nested pairings is counted as hairpin loop. Structural types and 3′-end dimer tags are derived from these fixed rules and no folding predictions are made.
- **API:** `dnakit.secondary_structure.analyze_dot_bracket(strands[required], dot_bracket[required], three_prime_window[optional])`
- **Input:** Required DNA strand collection and dot-bracket of consistent length; multiple strands are separated by `+`, `() [] {} <>` is supported.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import analyze_dot_bracket

result = analyze_dot_bracket(
    (DNASequence("ATCCTAGTTATAGGAT"),),
    "((((((....))))))",
)
print(result.structure_type, result.base_pair_count)
print(result.stem_lengths, result.hairpin_loop_lengths)
```

- **Example results:**

```text
hairpin 6
(6,) (4,)
```

## 2) Pairing probability and window accessibility

- **Function:** Verify the dimension, symmetry and probability range of the pairing probability matrix, calculate the paired and unpaired probabilities of each position, and provide reliable input for the ensemble indicator.
- **Calculation method:** It is required that the square matrix is ​​symmetrical, the elements are at `[0,1]` and the sum of each row and the diagonal is 1; the diagonal `P(i,i)` is used as the unpaired probability, and the row sum minus the diagonal is used as the total paired probability. Window accessibility is the arithmetic average of the unpaired edge probabilities of each site within the window.
- **API:** `dnakit.secondary_structure.pair_probability_metrics(strands[required], probability_matrix[required], accessibility_window_size[optional])`
- **Input:** Required DNA strand set and symmetric square matrix; diagonal is unpaired probability, the sum after each row including diagonal must be 1; optional accessibility window length.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import pair_probability_metrics

matrix = (
    (0.7, 0.0, 0.0, 0.3),
    (0.0, 0.8, 0.2, 0.0),
    (0.0, 0.2, 0.8, 0.0),
    (0.3, 0.0, 0.0, 0.7),
)
result = pair_probability_metrics(
    (DNASequence("ACGT"),),
    matrix,
    accessibility_window_size=2,
)
print(tuple(round(value, 1) for value in result.pairing_probabilities_by_base))
print(tuple(round(value, 1) for value in result.unpaired_probabilities_by_base))
print(result.most_accessible_window_start)
```

- **Example results:**

```text
(0.3, 0.2, 0.2, 0.3)
(0.7, 0.8, 0.8, 0.7)
1
```

## 3) Collection defects of target structure

- **Function:** Calculate ensemble defects based on the target pairing relationship and ensemble pairing probability, and estimate the expected number and proportion of bases that are not correctly paired according to the target structure.
- **Calculation method:** Determine the target state for each site: `P(i,i)` when unpaired, `P(i,j)` when paired with `j`; the normalized ensemble defect is `Σ[1 − P(target state)] / N`.
- **API:** `dnakit.secondary_structure.ensemble_defect_from_probabilities(target[required], probabilities[required])`
- **Input:** Required `SecondaryStructureSummary` Target and `PairProbabilityResult` describing the same DNA strand.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import (
    analyze_dot_bracket,
    ensemble_defect_from_probabilities,
    pair_probability_metrics,
)

matrix = (
    (0.7, 0.0, 0.0, 0.3),
    (0.0, 0.8, 0.2, 0.0),
    (0.0, 0.2, 0.8, 0.0),
    (0.3, 0.0, 0.0, 0.7),
)
strands = (DNASequence("ACGT"),)
target = analyze_dot_bracket(strands, "(())")
probabilities = pair_probability_metrics(strands, matrix, accessibility_window_size=2)
print(ensemble_defect_from_probabilities(target, probabilities))
```

- **Example results:**

```text
0.75
```

## 4) Target structure thermodynamic probability

- **Function:** Calculate the theoretical probability of the structure in the thermodynamic ensemble based on the free energy and partition function of the target structure, and use it to compare the relative proportions of candidate structures.
- **Calculation method:** Calculates according to the Boltzmann relation `P(target) = exp[−(Gtarget − Gensemble)/(RT)]` and limits the floating point result to `[0,1]`; the function does not calculate the two free energies by itself.
- **API:** `dnakit.secondary_structure.target_structure_probability(target_free_energy_kcal_per_mol[required], ensemble_free_energy_kcal_per_mol[required], temperature_celsius[optional])`
- **Input:** Required target structure and kcal/mol free energy of the ensemble; optional 0–100 °C temperature.
- **Sample code:**

```python
from dnakit.secondary_structure import target_structure_probability

same_energy = target_structure_probability(-2.0, -2.0)
higher_target_energy = target_structure_probability(-1.0, -2.0)
print(same_energy, round(higher_target_energy, 6))
```

- **Example results:**

```text
1.0 0.197404
```

## 5) NUPACK passive availability check

- **Function:** Check whether the current environment can import and call the compatible NUPACK backend, return version and capability information, so that the structure calculation can clearly report availability before running.
- **Calculation method:** Only Python package location and package metadata reading are used to check the `nupack` module, version and installation location; no modules are imported, no scientific calculations are performed, and no numerical accuracy is verified.
- **API:** `dnakit.secondary_structure.probe_nupack()`
- **Input:** None.
- **Sample code:**

```python
from dnakit.secondary_structure import probe_nupack

status = probe_nupack()
print(status.available, status.version, status.metadata["import_executed"])
```

- **Current local example results:**

```text
False None False
```

## 6) NUPACK single complex set analysis

- **Role:** Call NUPACK to calculate the MFE structure, free energy, pairing probability or structure sampling for one or more DNA strands and convert it into a unified result object for DNAKit.
- **Calculation method:** Use temperature, `[Na⁺]+[K⁺]`, `Mg²⁺` and ensemble to construct the NUPACK DNA model, and then explicitly call the external NUPACK's `pfunc`, `mfe`, `pairs`, `subopt`, `sample`, structure probability and defect interfaces; DNAKit is only responsible for input restrictions, result analysis and source records.
- **API:** `dnakit.secondary_structure.NupackAdapter.analyze_complex(strands[required], conditions[optional], ensemble[optional], target_structure[optional], suboptimal_energy_gap_kcal_per_mol[optional], num_samples[optional], accessibility_window_size[optional])`
- **Input:** Required DNA strand collection; optional NUPACK condition, ensemble, target dot-bracket, 0–20 kcal/mol suboptimal difference, and 0–100000 samples.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import NupackAdapter, probe_nupack

status = probe_nupack()
if status.available:
    result = NupackAdapter().analyze_complex(
        (DNASequence("ATCCTAGTTATAGGAT"),),
        target_structure="((((((....))))))",
        num_samples=100,
    )
    print(result.mfe_structures[0].summary.dot_bracket)
    print(result.ensemble_free_energy_kcal_per_mol)
else:
    print("NUPACK unavailable")
```

- **Current local example results:**

```text
NUPACK unavailable
```

## 7) NUPACK tube multi-complex balance

- **Function:** Call NUPACK to analyze the equilibrium composition of multiple complexes in the tube at a given chain concentration and temperature, and return the concentration and mass conservation information of each complex.
- **Calculation method:** Pass the naming chain, feed concentration, maximum complex size and target complex to the external NUPACK `tube_analysis`, which will calculate the equilibrium concentration and pairing probability of each complex; DNAKit then calculates the target/non-target ratio based on the sum of the complex concentrations enumerated this time.
- **API:** `dnakit.secondary_structure.NupackAdapter.analyze_tube(strands[required], concentrations_molar[required], target_strand_names[required], conditions[optional], max_complex_size[optional])`
- **Input:** Required Named DNA strand map, bond set identical molarity map, and target complex strand name; maximum 20 strands, `max_complex_size` is 1–4.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.secondary_structure import NupackAdapter, probe_nupack

status = probe_nupack()
if status.available:
    result = NupackAdapter().analyze_tube(
        {"a": DNASequence("CCC"), "b": DNASequence("GGG")},
        {"a": 1e-6, "b": 1e-6},
        target_strand_names=("a", "b"),
        max_complex_size=2,
    )
    print(result.target_complex_concentration_molar)
    print(result.target_complex_fraction, result.non_target_complex_fraction)
else:
    print("NUPACK unavailable")
```

- **Current local example results:**

```text
NUPACK unavailable
```

The external function boundaries of NUPACK can be found in its [official analysis documentation](https://docs.nupack.org/analysis/), [utility function documentation](https://docs.nupack.org/utilities/) and [model documentation](https://docs.nupack.org/model/); its separate [download and licensing requirements](https://www.nupack.org/download/overview) must be followed before installation.

!!! warning
    `probe_nupack().available=False`, dot-bracket comments, or Primer3 results must not be written as NUPACK calculation results.
