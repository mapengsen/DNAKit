# Physical and chemical properties

Focuses on the physicochemical properties of DNA sequences, as well as double-strand thermodynamics, base stacking, stability, hairpins, and dimer formation tendencies.

For the paper sources, internal rules and applicable boundaries of each calculation, see [FAQ: Calculation basis and references for physical and chemical properties ](../../faq.md#physicochemical-references).

The internal model and the optional external backend have different domains; all condition concentration units are mol/L.

## 1) molecular weight {#1-thermo-001-molecular-weight}

- **Purpose:** Estimate the theoretical molecular weight of DNA based on base composition, single- and double-stranded types, and end settings, and return the Da value, which can be used for molar concentration, mass concentration, and sample dosage conversion.
- **Calculation method:** Sum the masses of the A/C/G/T anhydrous deoxynucleotide residues for each strand, subtract the unphosphorylated end correction `61.96 Da`; add `79.0 Da` if 5′ phosphorylated, and add the mass of the complete reverse complement for duplex mode.
- **API:** `dnakit.thermodynamics.molecular_weight(sequence[required], strand[optional], five_prime_phosphorylated[optional], max_sequence_length[optional])`
- **Input:** Required linear canonical `DNASequence`; optional single/double stranded and 5′ phosphorylation state.
- **Output:** Returns a `MolecularWeightResult`; the example shows `value_dalton`, `value_kilodalton`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import molecular_weight

result = molecular_weight(DNASequence("ACGT"))
print(result.value_dalton, result.value_kilodalton)
```

- **Example output:**

```text
1173.84 1.17384
```

## 2) 260 nm extinction coefficient {#2-thermo-014-260-nm-extinction-coefficient}

- **Purpose:** Calculate the theoretical molar extinction coefficient `ε260` at 260 nm based on the single-stranded DNA sequence, which is used to convert nucleic acid concentration and absorbance according to Beer–Lambert's law.
- **Calculation method:** When the length is at least 2, sum the public extinction coefficients of all adjacent dinucleotides, and then subtract all internal single-base coefficients; single-nucleotides directly use the corresponding single-base coefficients.
- **API:** `dnakit.thermodynamics.extinction_coefficient_260nm(sequence[required], max_sequence_length[optional])`
- **Input:** Required 1–1,000,000 nt linear, Gap-free, A/C/G/T, single-stranded, unmodified `DNASequence`.
- **Output:** Returns an `ExtinctionCoefficientResult`; the example shows `value_m_inverse_cm_inverse`, `wavelength_nm`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import extinction_coefficient_260nm

result = extinction_coefficient_260nm(DNASequence("ACGT"))
print(result.value_m_inverse_cm_inverse, result.wavelength_nm)
```

- **Example output:**

```text
40300.0 260
```

- **Based on:** [Oligonucleotide quantification instructions for IDT](https://sg.idtdna.com/page/support-and-education/decoded-plus/oligo-quantification-getting-it-right), [ACGT public examples and parameter tables](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/pcr/quantitation-of-oligos), Warshaw–Tinoco 1966 (DOI [`10.1016/0022-2836(66)90115-X`](https://doi.org/10.1016/0022-2836(66)90115-X)) and Cantor–Warshaw–Shapiro 1970 (DOI [`10.1002/bip.1970.360090909`](https://doi.org/10.1002/bip.1970.360090909)).

## 3) melting temperature Tm {#3-thermo-002-melting-temperature-tm}

- **Purpose:** Calculate the DNA melting temperature Tm based on sequence length, composition or nearest neighbor model and experimental conditions, which is used to compare double-strand stability and design oligonucleotide experimental conditions.
- **Calculation method:** `wallace` method uses `2 × (A+T) + 4 × (G+C)`; `nearest_neighbor` method uses SantaLucia 1998 adjacent base stacking, termini, symmetry, salt concentration and chain concentration parameters to calculate Tm.
- **API:** `dnakit.thermodynamics.melting_temperature(sequence[required], method[optional], conditions[optional], config[optional])`
- **Input:** Required canonical `DNASequence`; optional `method`, `ThermodynamicConditions` and NN configurations.
- **Output:** Returns a `MeltingTemperatureResult`; the example shows `tm_celsius`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import melting_temperature

result = melting_temperature(DNASequence("AACG"), method="wallace")
print(result.tm_celsius)
```

- **Example output:**

```text
12.0
```

## 4) salt concentration correction {#4-thermo-003-salt-concentration-correction}

- **Purpose:** Make empirical corrections to Tm based on the concentration of monovalent ions such as Na⁺, K⁺, etc., and return the temperature changes in the specified salt environment to facilitate comparison of different buffer conditions.
- **Calculation method:** Use SantaLucia 1998 monovalent salt entropy modification formula `ΔS_salt = 0.368 × (N - 1) × ln([Na⁺] + [K⁺])`, where `N` is the sequence length and the concentration unit is mol/L.
- **API:** `dnakit.thermodynamics.salt_correction(sequence_length[required], conditions[optional])`
- **Input:** Required sequence length; optional `ThermodynamicConditions` containing Na⁺, K⁺ concentrations.
- **Output:** Returns a `SaltCorrectionResult`; the example shows `delta_s_cal_per_k_mol`.
- **Example code:**

```python
from dnakit.thermodynamics import ThermodynamicConditions, salt_correction

conditions = ThermodynamicConditions(sodium_molar=0.05)
result = salt_correction(10, conditions=conditions)
print(round(result.delta_s_cal_per_k_mol, 3))
```

- **Example output:**

```text
-9.922
```

## 5) local melting characteristics {#5-thermo-012-local-melting-characteristics}

- **Purpose:** Calculate local Tm along the sequence sliding window, return the position, temperature and overall range of each window, used to find local areas with abnormal stability.
- **Calculation method:** Generate fixed windows by `window_size` and `step`, repeatedly calling the selected Wallace or SantaLucia nearest-neighbor Tm method for each window, and summarizing local minima, maxima, and positions.
- **API:** `dnakit.thermodynamics.window_tm(sequence[required], window_size[required], step[optional], method[optional], conditions[optional], config[optional], max_windows[optional])`
- **Input:** Required sequence and window size; optional step size, method, conditions, NN configuration and window upper limit.
- **Output:** Returns a `WindowTmResult`; the example shows `windows`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import window_tm

result = window_tm(DNASequence("AACGTT"), 4, step=2, method="wallace")
print([(item.start, item.end, item.tm_celsius) for item in result.windows])
```

- **Example output:**

```text
[(0, 4, 12.0), (2, 6, 12.0)]
```

## 6) Synthetic risk {#6-eval-013-synthetic-risk}

- **Purpose:** Screen sequences one by one according to GC extreme, same-base continuous, repeated and known rules, and return hit positions, risk items and transparent scores for rule-based pre-checking before synthesis.
- **Calculation method:** Calculate the five risk components of global GC, local GC abnormal window, longest same-base sequence, tandem duplication and inverted duplication respectively, truncate to `[0, 1]` and average with equal weights; scores `<0.2`, `0.2–<0.5`, `≥0.5` are marked as low, medium and high risk in turn.
- **API:** `dnakit.evaluation.evaluate_synthesis_risk(value[required], config[optional])`; `config` uses `dnakit.evaluation.SynthesisRiskConfig`.
- **Input:** Linear, DNA with no explicit gaps; optional thresholds and resource caps.
- **Output:** Returns `EvaluationReport` with aggregate `metrics` and per-record `entries`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.evaluation import evaluate_synthesis_risk

report = evaluate_synthesis_risk(DNASequence("G" * 40 + "AT" * 10))
entry = report.entries[0]
print(entry.metrics["risk_level"])
print(entry.metrics["risk_score"])
```

- **Example output:**

```text
medium
0.38666666666666666
```

<span id="thermodynamic-properties"></span>

## 7) thermodynamic parameters {#1-thermo-004-thermodynamic-parameters}

- **Purpose:** Calculate ΔH, ΔS, ΔG and Tm upon fully complementary duplex formation using nearest neighbor parameters and explicit experimental conditions for quantitative comparison of the thermodynamic stabilities of candidate duplexes.
- **Calculation method:** Sum up the adjacent stacking, terminal onset, symmetry and monovalent salt contributions of SantaLucia 1998 to obtain ΔH and ΔS respectively, and then use `ΔG = ΔH - TΔS` and the concentration-corrected thermodynamic equation to calculate Tm.
- **API:** `dnakit.thermodynamics.nearest_neighbor(sequence[required], complement[optional], conditions[optional], config[optional])`
- **Input:** A canonical sequence is required; fully complementary strands, conditions and NN parameter sets are optional.
- **Output:** Returns a `NearestNeighborResult`; the example shows `delta_h_kcal_per_mol`, `delta_s_cal_per_k_mol`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, nearest_neighbor

conditions = ThermodynamicConditions(sodium_molar=1.0, strand_concentration_molar=1e-6)
result = nearest_neighbor(DNASequence("GTGCAT"), conditions=conditions)
print(result.delta_h_kcal_per_mol, result.delta_s_cal_per_k_mol)
```

- **Example output:**

```text
-40.0 -111.3
```

## 8) Nearest-neighbor {#2-thermo-005-nearest-neighbor}

- **Purpose:** Split the double strand into adjacent base pair steps, calculate the contribution of each step to ΔH, ΔS and ΔG respectively, and summarize the entire sequence to facilitate locating stable or unstable fragments.
- **Calculation method:** Split the sequence into all adjacent dinucleotide steps, query the SantaLucia 1998 ΔH/ΔS parameter table item by item, and then add terminal, symmetry and salt corrections; this item uses the same `nearest_neighbor()` result as thermodynamic parameters, but emphasizes the step-by-step details.
- **API:** `dnakit.thermodynamics.nearest_neighbor(sequence[required], complement[optional], conditions[optional], config[optional])`
- **Input:** Required 2–60 nt canonical linear sequence; optional fully complementary strand, condition, and `NearestNeighborConfig`.
- **Output:** Returns a `NearestNeighborResult`; the example shows `stacking_steps`, `tm_celsius`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, nearest_neighbor

result = nearest_neighbor(
    DNASequence("GTGCAT"),
    conditions=ThermodynamicConditions(sodium_molar=1.0, strand_concentration_molar=1e-6),
)
print([step.top_5to3 for step in result.stacking_steps], round(result.tm_celsius, 3))
```

- **Example output:**

```text
['GT', 'TG', 'GC', 'CA', 'AT'] 9.517
```

## 9) Duplex stability {#3-thermo-006-duplex-stability}

- **Purpose:** Check whether two DNAs satisfy the complete Watson–Crick complementary relationship; calculate the double-stranded thermodynamic results after checking to verify the paired sequence and its stability.
- **Calculation method:** The `native` backend requires the two strands to be complete reverse complements and uses SantaLucia nearest-neighbor results, with `Tm > configured temperature` judged stable; when `primer3-cli` is explicitly selected, Primer3 `ntthal` estimates a heterodimer structure that allows mismatches and dangling ends.
- **API:** `dnakit.thermodynamics.duplex_stability(sequence_a[required], sequence_b[required], conditions[optional], config[optional], backend[optional], adapter[optional], max_loop[optional], output_structure[optional])`
- **Input:** Required sequences A, B; optional conditions, `native`/`primer3-cli` backend, adapter, max loop and structure output.
- **Output:** Returns a `DuplexStabilityResult`; the example shows `fully_complementary`, `stable_at_temperature`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import duplex_stability

result = duplex_stability(DNASequence("GTGCAT"), DNASequence("ATGCAC"))
print(result.fully_complementary, result.stable_at_temperature)
```

- **Example output:**

```text
True False
```

## 10) base stacking {#4-thermo-007-base-stacking}

- **Purpose:** List the enthalpy, entropy and free energy contributions of each adjacent base stacking step separately to explain which local sequences the thermodynamic results of the entire double strand come from.
- **Calculation method:** Query the ΔH and ΔS parameters of SantaLucia 1998 for each adjacent dinucleotide and calculate `ΔG = ΔH - TΔS` at the specified temperature and then sum; this interface does not add end, symmetry and salt corrections.
- **API:** `dnakit.thermodynamics.stacking_interactions(sequence[required], temperature_celsius[optional], config[optional])`
- **Input:** Required canonical `DNASequence`; optional temperature and NN configuration.
- **Output:** Returns a `StackingResult`; the example shows `steps`, `total_delta_h_kcal_per_mol`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import stacking_interactions

result = stacking_interactions(DNASequence("AA"))
print(result.steps[0].top_5to3, result.total_delta_h_kcal_per_mol)
```

- **Example output:**

```text
AA -7.9
```

## 11) Hairpin {#5-thermo-008-hairpin}

- **Purpose:** Call the explicit structure backend to predict whether a single sequence forms a hairpin, and return indicators such as whether the structure exists, Tm and ΔG, for screening self-folding that may affect the use of primers or oligonucleotides.
- **Calculation method:** Call the user-installed Primer3 `ntthal` in `HAIRPIN` mode to search for single-molecule hairpin structures and resolve their Tm, ΔG, ΔH and ΔS under given salt concentration, chain concentration, temperature and `max_loop` conditions.
- **API:** `dnakit.thermodynamics.probe_primer3(ntthal_path[optional], thermodynamic_parameters_path[optional])`, `dnakit.thermodynamics.Primer3CLIAdapter.hairpin(sequence[required], conditions[optional], max_loop[optional], output_structure[optional])`
- **Input:** Required 1–60 nt canonical sequence; optional conditions, `max_loop=1..30` and structure output.
- **Output:** `probe_primer3` returns backend availability; hairpin analysis reports structure formation and Tm.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import Primer3CLIAdapter

adapter = Primer3CLIAdapter(
    ntthal_path="/opt/primer3/src/ntthal",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
result = adapter.hairpin(DNASequence("CCCCCATCCGATCAGGGGG"))
print(result.structure_found, result.tm_celsius)
```

- **Example output:**

```text
Results vary with the user-installed Primer3 version, parameter files, and conditions.
```

## 12) Self-dimer {#6-thermo-009-self-dimer}

- **Purpose:** Evaluate whether two molecules of the same sequence may form a self-dimer and return the predicted structure, Tm and ΔG for screening the risk of primer self-pairing.
- **Calculation method:** Use the same sequence as two inputs to call the Primer3 `ntthal` `ANY` mode installed by the user, search for the self-dimer thermodynamic structure and analyze the relevant parameters under the specified experimental conditions.
- **API:** `dnakit.thermodynamics.Primer3CLIAdapter.self_dimer(sequence[required], conditions[optional], max_loop[optional], output_structure[optional])`
- **Input:** Required 1–60 nt canonical sequence; optional conditions, max loop, and structure output.
- **Output:** Returns a `Primer3ThermodynamicResult`; the example shows `structure_found`, `delta_g_kcal_per_mol`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import Primer3CLIAdapter

adapter = Primer3CLIAdapter(
    ntthal_path="/opt/primer3/src/ntthal",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
result = adapter.self_dimer(DNASequence("CCCCCATCCGATCAGGGGG"))
print(result.structure_found, result.delta_g_kcal_per_mol)
```

- **Example output:**

```text
Results vary with the user-installed Primer3 version, parameter files, and conditions.
```

## 13) Heterodimer {#7-thermo-010-heterodimer}

- **Purpose:** Evaluate whether heterodimers are likely to form between two different sequences, returning the predicted structure, Tm and ΔG for checking for untargeted pairings between primer pairs or oligonucleotides.
- **Calculation method:** Input two different sequences into the user-installed Primer3 `ntthal` `ANY` mode, search for the heterodimer thermodynamic structure and resolve the relevant parameters under specified experimental conditions.
- **API:** `dnakit.thermodynamics.Primer3CLIAdapter.heterodimer(sequence_a[required], sequence_b[required], conditions[optional], max_loop[optional], output_structure[optional])`
- **Input:** Two 1–60 nt canonical sequences are required; optional conditions, max loop, and structure output.
- **Output:** Returns a `Primer3ThermodynamicResult`; the example shows `structure_found`, `delta_g_kcal_per_mol`.
- **Example code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import Primer3CLIAdapter

adapter = Primer3CLIAdapter(
    ntthal_path="/opt/primer3/src/ntthal",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
result = adapter.heterodimer(DNASequence("GTGCAT"), DNASequence("ATGCAC"))
print(result.structure_found, result.delta_g_kcal_per_mol)
```

- **Example output:**

```text
Results vary with the user-installed Primer3 version, parameter files, and conditions.
```
