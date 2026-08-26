# 光学与浓度换算性质

根据 DNA 序列和吸光度计算光学性质、摩尔浓度、质量浓度、物质的量和质量，并支持修饰基团校正。

所有消光系数单位为 `M⁻¹·cm⁻¹`，浓度和体积分别使用 `mol/L` 与 `L`。

## 1) 260 nm 单链摩尔消光系数

- **作用：** 按序列和相邻碱基模型计算单链 DNA 在 260 nm 的理论摩尔消光系数，输出 `M⁻¹·cm⁻¹` 数值，供吸光度与浓度换算使用。
- **API：** `dnakit.thermodynamics.extinction_coefficient_260nm(sequence[必须], max_sequence_length[可选])`
- **输入：** 必填线性、单链、canonical `A/C/G/T` 序列；可选序列长度上限。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import extinction_coefficient_260nm

result = extinction_coefficient_260nm(DNASequence("ACGT"))
print(result.value_m_inverse_cm_inverse, result.wavelength_nm, result.method)
```

- **示例结果：**

```text
40300.0 260 nearest-neighbor-hypochromicity
```

- **限制：** 结果是未修饰单链 DNA 的理论值，不适用于双链、Gap、IUPAC 模糊碱基或修饰碱基。双链应使用 `optical_properties()`。

## 2) 单链/双链理论光学性质

- **作用：** 在一次调用中计算单链或双链 DNA 的 `ε260`、分子量及相关光学属性，为后续 OD260、质量和物质的量换算提供统一参数。
- **API：** `dnakit.thermodynamics.optical_properties(sequence[必须], strand_type[可选], complement[可选], duplex_method[可选], hypochromicity_fraction[可选], modifications[可选])`
- **输入：** 必填 canonical DNA 序列；可选单/双链、显式互补链、双链模型、hypochromicity 分数和修饰项。
- **示例代码：**

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

- **示例结果：**

```text
40300.0 1173.84
52800.0 average-dsdna-base-pair-extinction
68510.0 sequence-specific-strand-sum-with-explicit-hypochromicity
```

- **限制：** `average-base-pair` 使用每碱基对 `13200 M⁻¹·cm⁻¹` 的平均换算，不是序列特异的双链 hypochromicity 预测；`strand-sum-hypochromicity` 必须由调用方提供实验或文献分数。

## 3) 1 OD260 对应的 nmol 和质量

- **作用：** 根据 DNA 的消光系数和分子量，把 1 OD260 对应的核酸量换算为 nmol 与 µg，避免使用不区分序列的固定经验系数。
- **API：** `dnakit.thermodynamics.optical_properties(sequence[必须], strand_type[可选], complement[可选], duplex_method[可选], hypochromicity_fraction[可选], modifications[可选])`
- **输入：** 必填 DNA 序列；可选链类型、双链模型和显式修饰。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import optical_properties

result = optical_properties(DNASequence("ACGT"))
print(round(result.one_od260_nmol, 4), round(result.one_od260_microgram, 4))
```

- **示例结果：**

```text
24.8139 29.1275
```

- **限制：** 这里 1 OD260 定义为 1 cm 参考光程下 `A260 × 体积(mL) = 1`；其准确性直接取决于所选 ε260 模型和修饰参数。

## 4) A260 到摩尔浓度、质量浓度和总量

- **作用：** 使用 A260、光程、消光系数和样品体积计算摩尔浓度、质量浓度、总物质的量及总质量，用于定量核酸样品。
- **API：** `dnakit.thermodynamics.concentration_from_a260(measured_a260[必须], properties[必须], path_length_cm[可选], dilution_factor[可选], label_corrections[可选], volume_liter[可选])`
- **输入：** 必填实测 A260 和 `OpticalPropertiesResult`；可选光程、稀释倍数、染料校正和样品体积。
- **示例代码：**

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

- **示例结果：**

```text
0.403
10.0 11.7384
10.0 11.7384
```

- **限制：** 光程和稀释倍数必须显式且为正数；不提供 `volume_liter` 时，总物质的量和总质量为 `None`。

## 5) 浓度、物质的量和质量互换

- **作用：** 结合分子量和体积，在摩尔浓度、质量浓度、物质的量和质量之间进行一致换算，并明确返回各结果单位。
- **API：** `dnakit.thermodynamics.convert_oligo_quantity(molecular_weight_dalton[必须], volume_liter[可选], molar_concentration_molar[可选], mass_concentration_g_per_l[可选], amount_mol[可选], mass_g[可选])`
- **输入：** 必填分子量；四种浓度/总量输入中必须且只能提供一种，浓度转总量时还必须提供体积。
- **示例代码：**

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

- **示例结果：**

```text
10.0 11.7384
molar_concentration_molar
```

- **限制：** 该 API 只做分子量、体积和数量之间的换算，不会从序列自动推断分子量，也不校正测量误差。

## 6) 染料和修饰基团的显式校正

- **作用：** 根据用户提供的染料、末端基团或其他修饰增量，校正 `ε260`、分子量和 A260 派生结果，使带修饰寡核苷酸的浓度换算可审计。
- **API：** `dnakit.thermodynamics.OpticalModification(name[必须], count[可选], extinction_coefficient_260_delta_m_inverse_cm_inverse[可选], molecular_weight_delta_dalton[可选])`、`dnakit.thermodynamics.LabelAbsorbanceCorrection(name[必须], absorbance_at_label_max[必须], a260_correction_factor[必须])`
- **输入：** 修饰名称及显式的修正值；`OpticalModification` 传入 `optical_properties()`，`LabelAbsorbanceCorrection` 传入 `concentration_from_a260()`。
- **示例代码：**

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

- **示例结果：**

```text
42300.0 1373.84
0.02 0.48 11.3475
```

- **限制：** DNAKit 不内置厂商或批次特异的染料参数；Δε、ΔMW 和交叉吸收因子必须由调用方从实验或可追溯资料中提供。
