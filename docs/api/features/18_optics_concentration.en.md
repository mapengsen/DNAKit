# Conversion

Calculate optical properties, molar concentration, mass concentration, amount and mass of substances based on DNA sequence and absorbance, and support modification group correction.

All extinction coefficient units are `M⁻¹·cm⁻¹`, concentration and volume use `mol/L` and `L` respectively.

For the paper sources and internal formulas of various calculations, see [FAQ: Calculation basis and references for optical and concentration conversions ](../../faq.md#optics-concentration-references).

## 1) Single chain molar extinction coefficient at 260 nm

- **Function:** Calculate the theoretical molar extinction coefficient of single-stranded DNA at 260 nm based on the sequence and adjacent base models, and output the `M⁻¹·cm⁻¹` value for use in absorbance and concentration conversion.
- **Calculation method:** Look up the single-base `ε260` parameter for one nucleotide; for lengths of at least 2, calculate `ε260 = Σ adjacent-dinucleotide coefficients − Σ internal single-base coefficients` deterministically.
- **API:** `dnakit.thermodynamics.extinction_coefficient_260nm(sequence[required], max_sequence_length[optional])`
- **Input:** Required linear, single-stranded, canonical `A/C/G/T` sequence; optional upper sequence length.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import extinction_coefficient_260nm

result = extinction_coefficient_260nm(DNASequence("ACGT"))
print(result.value_m_inverse_cm_inverse, result.wavelength_nm, result.method)
```

- **Example results:**

```text
40300.0 260 nearest-neighbor-hypochromicity
```

## 2) Single chain/double chain theoretical optical properties

- **Function:** Calculate `ε260`, molecular weight and related optical properties of single-stranded or double-stranded DNA in one call, providing unified parameters for subsequent OD260, mass and substance quantity conversion.
- **Calculation method:** Single-stranded calculations use adjacent-dinucleotide `ε260` and the molecular-weight formula for anhydrous DNA residues; double-stranded `average-base-pair` uses `13200 × base-pair count` as `ε260`, while `strand-sum-hypochromicity` uses `(ε1 + ε2) × (1 − h)`, where `h` is supplied by the caller. Modification `Δε` and `ΔMW` values are finally added by count.
- **API:** `dnakit.thermodynamics.optical_properties(sequence[required], strand_type[optional], complement[optional], duplex_method[optional], hypochromicity_fraction[optional], modifications[optional])`
- **Input:** Required canonical DNA sequence; optional single/double strand, explicit complementary strand, double strand model, hypochromicity score, and modifications.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import optical_properties

single = optical_properties(DNASequence("ACGT"))
average = optical_properties(DNASequence("ACGT"), strand_type="double")
explicit = optical_properties(
    DNASequence("ACGT"),
    strand_type="double",
    duplex_method="strand-sum-hypochromicity",
    hypochromicity_fraction=0.15,
)
print(
    single.extinction_coefficient_260_m_inverse_cm_inverse,
    round(single.molecular_weight_dalton, 2),
)
print(average.extinction_coefficient_260_m_inverse_cm_inverse, average.method)
print(explicit.extinction_coefficient_260_m_inverse_cm_inverse, explicit.method)
```

- **Example results:**

```text
40300.0 1173.84
52800.0 average-dsdna-base-pair-extinction
68510.0 sequence-specific-strand-sum-with-explicit-hypochromicity
```

## 3) nmol and mass corresponding to 1 OD260

- **Function:** Based on the extinction coefficient and molecular weight of DNA, convert the amount of nucleic acid corresponding to 1 OD260 into nmol and µg to avoid using fixed empirical coefficients that do not distinguish sequences.
- **Calculation method:** Defined by `1 OD260 = A260 × volume (mL) = 1` at a 1 cm reference optical path length, converted using `nmol/OD260 = 10⁶/ε260` and `µg/OD260 = 1000 × MW/ε260`.
- **API:** `dnakit.thermodynamics.optical_properties(sequence[required], strand_type[optional], complement[optional], duplex_method[optional], hypochromicity_fraction[optional], modifications[optional])`
- **Input:** Required DNA sequence; optional strand type, duplex model, and explicit modifications.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import optical_properties

result = optical_properties(DNASequence("ACGT"))
print(round(result.one_od260_nmol, 4), round(result.one_od260_microgram, 4))
```

- **Example results:**

```text
24.8139 29.1275
```

## 4) A260 to molar concentration, mass concentration and total amount

- **Function:** Use A260, optical path, extinction coefficient and sample volume to calculate molar concentration, mass concentration, amount of total substance and total mass for quantification of nucleic acid samples.
- **Calculation method:** First calculate `A260_corrected = A260_measured − Σ(label peak × cross-absorption factor)`, then obtain molar concentration from the Beer–Lambert law as `c = A260_corrected × dilution factor / (ε260 × path length)`; mass concentration is `c × MW`, and when volume is provided, calculate `n = cV` and total mass.
- **API:** `dnakit.thermodynamics.concentration_from_a260(measured_a260[required], properties[required], path_length_cm[optional], dilution_factor[optional], label_corrections[optional], volume_liter[optional])`
- **Input:** Required measured A260 and `OpticalPropertiesResult`; optional path length, dilution factor, dye calibration and sample volume.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import concentration_from_a260, optical_properties

properties = optical_properties(DNASequence("ACGT"))
result = concentration_from_a260(
    0.403,
    properties,
    path_length_cm=1.0,
    dilution_factor=1.0,
    volume_liter=0.001,
)
assert result.amount_mol is not None and result.mass_microgram is not None
print(result.corrected_a260)
print(
    round(result.molar_concentration_micromolar, 4),
    round(result.mass_concentration_ng_per_microliter, 4),
)
print(round(result.amount_mol * 1e9, 4), round(result.mass_microgram, 4))
```

- **Example results:**

```text
0.403
10.0 11.7384
10.0 11.7384
```

## 5) Concentration, quantity and mass interchange of substances

- **What it does:** Combines molecular weight and volume to perform consistent conversions between molar concentration, mass concentration, amount and mass of a substance, and returns each result unit explicitly.
- **Calculation method:** Use `m = n × MW`, `n = m/MW`, `c = n/V` and `n = cV` for pure dimensional conversion; mass concentration and molar concentration are interchanged through `MW`, and no model estimation is performed.
- **API:** `dnakit.thermodynamics.convert_oligo_quantity(molecular_weight_dalton[required], volume_liter[optional], molar_concentration_molar[optional], mass_concentration_g_per_l[optional], amount_mol[optional], mass_g[optional])`
- **Input:** Required molecular weight; only one of the four concentration/total inputs must be provided, and the volume must also be provided when converting concentration to total amount.
- **Sample code:**

```python
from dnakit.thermodynamics import convert_oligo_quantity

result = convert_oligo_quantity(
    1173.84,
    volume_liter=0.001,
    molar_concentration_molar=1e-5,
)
print(round(result.amount_nmol, 4), round(result.mass_microgram, 4))
print(result.input_kind)
```

- **Example results:**

```text
10.0 11.7384
molar_concentration_molar
```

## 6) Explicit correction of dyes and modifying groups

- **Function:** Correct `ε260`, molecular weight, and A260 derived results based on user-supplied dye, end group, or other modification increments, making the concentration conversion of modified oligonucleotides auditable.
- **Calculation method:** For each modification, sum `count × Δε260` and `count × ΔMW` and add them to the unmodified value; dye cross-absorption is subtracted from measured A260 as `label peak × correction factor`. All correction parameters are supplied by the caller.
- **API:** `dnakit.thermodynamics.OpticalModification(name[required], count[optional], extinction_coefficient_260_delta_m_inverse_cm_inverse[optional], molecular_weight_delta_dalton[optional])`, `dnakit.thermodynamics.LabelAbsorbanceCorrection(name[required], absorbance_at_label_max[required], a260_correction_factor[required])`
- **Input:** Modifier name and explicit modifier value; `OpticalModification` passes in `optical_properties()`, `LabelAbsorbanceCorrection` passes in `concentration_from_a260()`.
- **Sample code:**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import (
    LabelAbsorbanceCorrection,
    OpticalModification,
    concentration_from_a260,
    optical_properties,
)

properties = optical_properties(
    DNASequence("ACGT"),
    modifications=(
        OpticalModification(
            "fluorophore",
            count=2,
            extinction_coefficient_260_delta_m_inverse_cm_inverse=1000.0,
            molecular_weight_delta_dalton=100.0,
        ),
    ),
)
result = concentration_from_a260(
    0.5,
    properties,
    label_corrections=(
        LabelAbsorbanceCorrection(
            "fluorophore",
            absorbance_at_label_max=0.2,
            a260_correction_factor=0.1,
        ),
    ),
)
print(
    properties.extinction_coefficient_260_m_inverse_cm_inverse,
    round(properties.molecular_weight_dalton, 2),
)
print(
    round(result.label_a260_subtracted, 4),
    result.corrected_a260,
    round(result.molar_concentration_micromolar, 4),
)
```

- **Example results:**

```text
42300.0 1373.84
0.02 0.48 11.3475
```
