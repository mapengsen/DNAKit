# Double chain thermodynamics extension

Calculate thermodynamics, stability, binding equilibria, melting curves, end stability, and cosolvent corrections for fully complementary DNA duplexes.

ΔH, ΔS, specified temperature ΔG and Tm for duplex formation using versioned SantaLucia 1998 nearest neighbor parameters.

For the paper sources and internal formulas of various calculations, see [FAQ: Calculation basis and references for double-chain thermodynamic expansion ](../../faq.md#duplex-thermodynamics-references).

## 1) ΔH, ΔS, ΔG and Tm of completely complementary double strands

- **Function:** Calculate ΔH, ΔS, ΔG and Tm of a fully complementary DNA duplex using the nearest neighbor model to quantify duplex stability at given salt concentration, strand concentration and temperature.
- **Calculation method:** Accumulate adjacent stacking, initiation at both ends and self-complementary symmetry terms according to SantaLucia 1998 parameters, and add `0.368 × (N−1) × ln([Na⁺]+[K⁺])` to `ΔS`; then use `ΔG = ΔH − TΔS/1000` and `Tm = 1000ΔH/(ΔS + R ln(Ct/divisor)) − 273.15`, where self-complementary `divisor=1`, otherwise 4.
- **API:** `dnakit.thermodynamics.nearest_neighbor(sequence[required], complement[optional], conditions[optional], config[optional])`
- **Input:** Required 2–60 nt canonical linear DNA; optional full reverse complementary strand, temperature, Na⁺/K⁺ total salt, strand concentration, and NN parameter configuration.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, nearest_neighbor

conditions = ThermodynamicConditions(
    temperature_celsius=37.0,
    sodium_molar=1.0,
    strand_concentration_molar=1e-6,
)
result = nearest_neighbor(DNASequence("GTGCAT"), conditions=conditions)
print(result.delta_h_kcal_per_mol, result.delta_s_cal_per_k_mol)
print(round(result.delta_g_kcal_per_mol, 4), round(result.tm_celsius, 4))
```

- **Example results:**

```text
-40.0 -111.3
-5.4803 9.5175
```

## 2) Unify double-strand stability results

- **Function:** Check the complementary relationship between two sequences and summarize stability indicators such as ΔG and Tm to facilitate judging whether the candidate pairing is stable enough under specified conditions.
- **Calculation method:** The `native` path requires the second strand to be a complete reverse complement and directly summarizes the nearest-neighbor results above; the `primer3-cli` path calls the user-installed `ntthal` to select a mismatch/dangling-end structure. Stable booleans are determined by `Tm > configured temperature`; Primer3 paths also require a structure to be found.
- **API:** `dnakit.thermodynamics.duplex_stability(sequence_a[required], sequence_b[required], conditions[optional], config[optional], backend[optional], adapter[optional], max_loop[optional], output_structure[optional])`
- **Input:** Two DNA sequences are required; optional conditions, native configuration, or explicit Primer3 adapter.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, duplex_stability

result = duplex_stability(
    DNASequence("GTGCAT"),
    DNASequence("ATGCAC"),
    conditions=ThermodynamicConditions(
        sodium_molar=1.0,
        strand_concentration_molar=1e-6,
    ),
)
print(
    result.fully_complementary,
    round(result.delta_g_kcal_per_mol, 4),
    round(result.tm_celsius, 4),
)
```

- **Example results:**

```text
True -5.4803 9.5175
```

## 3) Contribution of adjacent base pairs to steps

- **Effect:** Step by step calculate the contribution of each adjacent base pair to ΔH, ΔS and ΔG, returning a decomposition of band positions to account for the source of local stability in the duplex.
- **Calculation method:** Look up SantaLucia `ΔH/ΔS` table for each adjacent two-base step in the sequence, calculate `ΔG = ΔH − TΔS/1000` step by step and sum; this decomposition does not add terminal, symmetry and salt correction terms.
- **API:** `dnakit.thermodynamics.stacking_interactions(sequence[required], temperature_celsius[optional], config[optional])`
- **Input:** Required canonical linear DNA of at least 2 nt; optional temperature and NN parameter configuration.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import stacking_interactions

result = stacking_interactions(DNASequence("GTGCAT"), temperature_celsius=37.0)
first = result.steps[0]
print(len(result.steps), first.top_5to3, first.bottom_3to5)
print(first.delta_h_kcal_per_mol, first.delta_s_cal_per_k_mol)
print(round(result.total_delta_g_kcal_per_mol, 4))
```

- **Example results:**

```text
5 GT CA
-8.4 -22.4
-7.4771
```

## 4) Conditions and Na⁺/K⁺ unit price salt

- **Function:** Use structured objects to uniformly record temperature, salt concentration, chain concentration and co-solvent conditions, so that different thermodynamic calculations can be reused and accurately compared in the same experimental environment.
- **Calculation method:** This object mainly performs value range checking and condition recording; the native nearest neighbor model only uses `[Na⁺]+[K⁺]` as the total unit price salt. The `Mg²⁺`, dNTP, DMSO and formamide fields are only passed to external backends or correction functions that explicitly support them.
- **API:** `dnakit.thermodynamics.ThermodynamicConditions(temperature_celsius[optional], sodium_molar[optional], potassium_molar[optional], magnesium_molar[optional], dntp_molar[optional], strand_concentration_molar[optional], dmso_percent[optional], dmso_factor_celsius_per_percent[optional], formamide_molar[optional], salt_model[optional])`
- **Input:** Use `mol/L` for all concentrations; °C for temperatures.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, nearest_neighbor

sodium = nearest_neighbor(
    DNASequence("GTGCAT"),
    conditions=ThermodynamicConditions(sodium_molar=0.05, potassium_molar=0.0),
)
potassium = nearest_neighbor(
    DNASequence("GTGCAT"),
    conditions=ThermodynamicConditions(sodium_molar=0.0, potassium_molar=0.05),
)
print(round(sodium.tm_celsius, 4), round(potassium.tm_celsius, 4))
print(sodium.conditions.monovalent_molar)
```

- **Example results:**

```text
-6.0845 -6.0845
0.05
```

## 5) Ka, Kd and double chain ratio

- **Function:** Calculate the equilibrium constants Ka, Kd and the expected double-chain ratio based on the standard free energy ΔG and chain concentration, which are used to convert thermodynamic stability into equilibrium binding amount.
- **Calculation method:** First obtain `ΔG` from the nearest neighbor model, and then calculate according to `Ka = exp(−ΔG/RT)` and `Kd = 1/Ka`. The double-strand concentration and ratio are obtained from the positive roots of the quadratic equation of the ideal two-state mass conservation equation, and the self-complementary and non-self-complementary cases are distinguished.
- **API:** `dnakit.thermodynamics.binding_equilibrium(sequence[required], complement[optional], conditions[optional], config[optional])`
- **Input:** Required canonical DNA; optional fully complementary strand, temperature, monovalent salt, total oligonucleotide strand concentration, and NN configuration.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, binding_equilibrium

result = binding_equilibrium(
    DNASequence("GTGCAT"),
    conditions=ThermodynamicConditions(
        sodium_molar=1.0,
        strand_concentration_molar=1e-6,
    ),
)
print(f"{result.association_constant_m_inverse:.6e}")
print(f"{result.dissociation_constant_molar:.6e}")
print(round(result.duplex_fraction, 6), result.self_complementary)
```

- **Example results:**

```text
7.272209e+03
1.375098e-04
0.00361 False
```

## 6) Theoretical melting curve

- **Function:** Calculate the double-strand ratio point by point in a given temperature range, and return the theoretical melting curve and transition area, which is used to observe the double-strand dissociation trend when the temperature increases.
- **Calculation method:** Recalculate the nearest neighbors `ΔG`, `Ka` and the ideal two-state double chain ratio for each input temperature; if two adjacent points span 0.5, linear interpolation is used to obtain the curve midpoint temperature.
- **API:** `dnakit.thermodynamics.theoretical_melting_curve(sequence[required], temperatures_celsius[required], complement[optional], conditions[optional], config[optional], progress[optional])`
- **Input:** Required DNA sequence and 2–100001 strictly increasing temperature points at 0–100 °C; optional conditions and progress callbacks.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, theoretical_melting_curve

progress = []
result = theoretical_melting_curve(
    DNASequence("GTGCAT"),
    range(0, 51, 5),
    conditions=ThermodynamicConditions(
        sodium_molar=1.0,
        strand_concentration_molar=1e-6,
    ),
    progress=lambda completed, total: progress.append((completed, total)),
)
print(
    round(result.points[0].duplex_fraction, 6),
    round(result.points[-1].duplex_fraction, 6),
)
print(round(result.midpoint_temperature_celsius or 0.0, 4), progress[-1])
```

- **Example results:**

```text
0.815339 0.000267
9.4772 (11, 11)
```

## 7) 5′/3′ end stability

- **Function:** Calculate the free energy or related stability index of the specified window at both ends of DNA respectively, and return the endpoint difference, which is used to analyze end stability asymmetry and primer 3′ end characteristics.
- **Calculation method:** Take the 5′ end and 3′ end windows of equal length, and call the fully complementary nearest neighbor model respectively; compare the two `ΔG`, and the end with higher value (less negative) is marked as less stable.
- **API:** `dnakit.thermodynamics.terminal_stability(sequence[required], window_size[optional], conditions[optional], config[optional])`
- **Input:** Required 2–60 nt canonical DNA; optional 2 to Window size, conditions and NN configuration within sequence length range.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, terminal_stability

result = terminal_stability(
    DNASequence("AACCGGTT"),
    window_size=5,
    conditions=ThermodynamicConditions(
        sodium_molar=1.0,
        strand_concentration_molar=1e-6,
    ),
)
print(result.five_prime_sequence, round(result.five_prime_delta_g_kcal_per_mol, 4))
print(result.three_prime_sequence, round(result.three_prime_delta_g_kcal_per_mol, 4))
print(result.less_stable_end)
```

- **Example results:**

```text
AACCG -4.4624
CGGTT -4.4624
equal
```

## 8) DMSO and formamide experience correction

- **Effect:** Apply an explicit empirical correction to the base Tm based on DMSO or formamide concentration, returning the correction value and correction amount for approximate comparison of cosolvent-containing conditions.
- **Calculation method:** Calculate `ΔTm_DMSO = −factor × DMSO%` and `ΔTm_formamide = (0.453 × GC fraction − 2.88) × formamide (mol/L)` according to the additive empirical formula in the Primer3 manual, and then add them to the uncorrected Tm.
- **API:** `dnakit.thermodynamics.cosolvent_tm_correction(sequence[required], uncorrected_tm_celsius[required], dmso_percent[optional], dmso_factor_celsius_per_percent[optional], formamide_molar[optional])`
- **Input:** Required canonical DNA and uncorrected Tm; optional DMSO volume percent, correction factor per percent, and formamide molar concentration.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import cosolvent_tm_correction

result = cosolvent_tm_correction(
    DNASequence("ACGT"),
    60.0,
    dmso_percent=5.0,
    formamide_molar=1.0,
)
print(result.dmso_delta_tm_celsius)
print(round(result.formamide_delta_tm_celsius, 4))
print(round(result.corrected_tm_celsius, 4))
```

- **Example results:**

```text
-3.0
-2.6535
54.3465
```

## 9) Mg²⁺, dNTP, mismatch and dangling end for Primer3 CLI

- **Use:** Compute Tm, hairpin, autodimer or heterodimer results under complex salt conditions through the explicit Primer3 backend to supplement structural thermodynamic analysis not covered by the native model.
- **Calculation method:** DNAKit only verifies the input, assembly parameters and parses the output; Tm is calculated by Primer3 `oligotm` explicitly configured by the user, and hairpin/dimer is calculated by `ntthal`. The specific thermodynamic model, Mg²⁺/dNTP handling, and mismatch/dangling-end structure choices depend on the actual Primer3 version and parameter catalog.
- **API:** `dnakit.thermodynamics.probe_primer3(oligotm_path[optional], ntthal_path[optional], thermodynamic_parameters_path[optional])`, `dnakit.thermodynamics.Primer3CLIAdapter.tm(sequence[required], conditions[optional])`, `dnakit.thermodynamics.Primer3CLIAdapter.heterodimer(sequence_a[required], sequence_b[required], conditions[optional], max_loop[optional], output_structure[optional])`
- **Input:** Tm required 2–36 nt canonical DNA; heterodimer required two 1–60 nt canonical DNA; optional Na⁺/K⁺, Mg²⁺, dNTP, strand concentration and structure options.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import (
    Primer3CLIAdapter,
    ThermodynamicConditions,
)

conditions = ThermodynamicConditions(
    sodium_molar=0.05,
    magnesium_molar=0.0015,
    dntp_molar=0.0006,
    strand_concentration_molar=50e-9,
)
adapter = Primer3CLIAdapter(
    oligotm_path="/opt/primer3/src/oligotm",
    ntthal_path="/opt/primer3/src/ntthal",
    thermodynamic_parameters_path="/opt/primer3/src/primer3_config",
)
tm = adapter.tm(DNASequence("GTAAAACGACGGCCAGT"), conditions=conditions)
heterodimer = adapter.heterodimer(
    DNASequence("GTGCAT"),
    DNASequence("ATGCAC"),
    conditions=conditions,
)
print(tm.tm_celsius, heterodimer.delta_g_kcal_per_mol)
```

- **Current local example results:**

```text
Results vary with the user-installed Primer3 version, parameter files, and conditions.
```

!!! warning "Arbitrarily modified base thermodynamics are not currently supported"
    Modifying bases requires a parameter set that matches the specific chemical structure and adjacent environment; there is currently no public call entry and will not be silently replaced with canonical parameters.
