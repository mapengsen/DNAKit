# 双链热力学扩展

计算完全互补 DNA 双链的热力学、稳定性、结合平衡、熔解曲线、末端稳定性和共溶剂修正。

双链形成的 ΔH、ΔS、指定温度 ΔG 和 Tm 使用版本化 SantaLucia 1998 最近邻参数。

## 1) 完全互补双链的 ΔH、ΔS、ΔG 和 Tm

- **作用：** 使用最近邻模型计算完全互补 DNA 双链的 ΔH、ΔS、ΔG 和 Tm，量化给定盐浓度、链浓度及温度下的双链稳定性。
- **API：** `dnakit.thermodynamics.nearest_neighbor(sequence[必须], complement[可选], conditions[可选], config[可选])`
- **输入：** 必填 2–60 nt canonical 线性 DNA；可选完全反向互补链、温度、Na⁺/K⁺ 总单价盐、链浓度和 NN 参数配置。
- **示例代码：**

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

- **示例结果：**

```text
-40.0 -111.3
-5.4803 9.5175
```

- **限制：** native 模型只支持线性、canonical、完全互补 DNA；不支持 mismatch、dangling end、Mg²⁺、dNTP 或修饰碱基参数。

## 2) 统一双链稳定性结果

- **作用：** 检查两条序列的互补关系并汇总 ΔG、Tm 等稳定性指标，便于判断候选配对在指定条件下是否足够稳定。
- **API：** `dnakit.thermodynamics.duplex_stability(sequence_a[必须], sequence_b[必须], conditions[可选], config[可选], backend[可选], adapter[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填两条 DNA 序列；可选条件、native 配置或显式 Primer3 adapter。
- **示例代码：**

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

- **示例结果：**

```text
True -5.4803 9.5175
```

- **限制：** 默认 `backend="native"` 要求完全反向互补；只有用户提供含显式 `ntthal_path` 的 adapter 并选择 `backend="primer3-cli"`，才会处理 Primer3 选中的 mismatch/dangling-end 结构。

## 3) 相邻碱基对步骤贡献

- **作用：** 逐步计算每一对相邻碱基对对 ΔH、ΔS 和 ΔG 的贡献，返回带位置的分解结果，用于解释双链中的局部稳定性来源。
- **API：** `dnakit.thermodynamics.stacking_interactions(sequence[必须], temperature_celsius[可选], config[可选])`
- **输入：** 必填至少 2 nt 的 canonical 线性 DNA；可选温度和 NN 参数配置。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import stacking_interactions

result = stacking_interactions(DNASequence("GTGCAT"), temperature_celsius=37.0)
first = result.steps[0]
print(len(result.steps), first.top_5to3, first.bottom_3to5)
print(first.delta_h_kcal_per_mol, first.delta_s_cal_per_k_mol)
print(round(result.total_delta_g_kcal_per_mol, 4))
```

- **示例结果：**

```text
5 GT CA
-8.4 -22.4
-7.4771
```

- **限制：** 步骤合计不包含 initiation、symmetry 和 salt 项，因此不等同于完整双链热力学总量。

## 4) 条件与 Na⁺/K⁺ 单价盐

- **作用：** 用结构化对象统一记录温度、盐浓度、链浓度和共溶剂条件，使不同热力学计算能够复用并准确比较同一实验环境。
- **API：** `dnakit.thermodynamics.ThermodynamicConditions(temperature_celsius[可选], sodium_molar[可选], potassium_molar[可选], magnesium_molar[可选], dntp_molar[可选], strand_concentration_molar[可选], dmso_percent[可选], dmso_factor_celsius_per_percent[可选], formamide_molar[可选], salt_model[可选])`
- **输入：** 所有浓度使用 `mol/L`；温度使用 °C。
- **示例代码：**

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

- **示例结果：**

```text
-6.0845 -6.0845
0.05
```

- **限制：** native NN 不会把 Mg²⁺/dNTP 静默折算为单价盐；这些字段只供显式支持它们的后端使用。

## 5) Ka、Kd 和双链比例

- **作用：** 根据标准自由能 ΔG 和链浓度计算平衡常数 Ka、Kd 及预计双链比例，用于把热力学稳定性转换为平衡结合量。
- **API：** `dnakit.thermodynamics.binding_equilibrium(sequence[必须], complement[可选], conditions[可选], config[可选])`
- **输入：** 必填 canonical DNA；可选完全互补链、温度、单价盐、总寡核苷酸链浓度和 NN 配置。
- **示例代码：**

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

- **示例结果：**

```text
7.272209e+03
1.375098e-04
0.00361 False
```

- **限制：** 这是理想两态、标准态平衡模型，不包含多中间态、竞争复合物、动力学或活度系数。

## 6) 理论熔解曲线

- **作用：** 在给定温度范围逐点计算双链比例，返回理论熔解曲线及转变区域，用于观察温度升高时双链解离趋势。
- **API：** `dnakit.thermodynamics.theoretical_melting_curve(sequence[必须], temperatures_celsius[必须], complement[可选], conditions[可选], config[可选], progress[可选])`
- **输入：** 必填 DNA 序列和 2–100001 个严格递增、位于 0–100 °C 的温度点；可选条件和进度回调。
- **示例代码：**

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

- **示例结果：**

```text
0.815339 0.000267
9.4772 (11, 11)
```

- **限制：** 曲线是理想平衡双链分数，不是仪器吸光度响应、动力学轨迹或热容模型。

## 7) 5′/3′ 末端稳定性

- **作用：** 分别计算 DNA 两端指定窗口的自由能或相关稳定性指标，返回端点差异，用于分析末端稳定性不对称和引物 3′ 端特征。
- **API：** `dnakit.thermodynamics.terminal_stability(sequence[必须], window_size[可选], conditions[可选], config[可选])`
- **输入：** 必填 2–60 nt canonical DNA；可选 2 至序列长度范围内的窗口大小、条件和 NN 配置。
- **示例代码：**

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

- **示例结果：**

```text
AACCG -4.4624
CGGTT -4.4624
equal
```

- **限制：** 它只比较两个等长窗口，不是整条序列的二级结构、末端 fraying 或反应动力学预测。

## 8) DMSO 和甲酰胺经验修正

- **作用：** 根据 DMSO 或甲酰胺浓度对基础 Tm 施加显式经验修正，返回校正值和修正量，用于近似比较含共溶剂条件。
- **API：** `dnakit.thermodynamics.cosolvent_tm_correction(sequence[必须], uncorrected_tm_celsius[必须], dmso_percent[可选], dmso_factor_celsius_per_percent[可选], formamide_molar[可选])`
- **输入：** 必填 canonical DNA 和未修正 Tm；可选 DMSO 体积百分比、每百分比校正因子和甲酰胺摩尔浓度。
- **示例代码：**

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

- **示例结果：**

```text
-3.0
-2.6535
54.3465
```

- **限制：** 这是 Primer3 手册形式的经验加和修正，不会重新计算含共溶剂的 ΔH/ΔS，也不是机理性自由能模型。

## 9) Primer3 CLI 的 Mg²⁺、dNTP、mismatch 和 dangling end

- **作用：** 通过显式 Primer3 后端计算复杂盐条件下的 Tm、发卡、自二聚体或异二聚体结果，用于补充原生模型未覆盖的结构型热力学分析。
- **API：** `dnakit.thermodynamics.probe_primer3(oligotm_path[可选], ntthal_path[可选], thermodynamic_parameters_path[可选])`、`dnakit.thermodynamics.Primer3CLIAdapter.tm(sequence[必须], conditions[可选])`、`dnakit.thermodynamics.Primer3CLIAdapter.heterodimer(sequence_a[必须], sequence_b[必须], conditions[可选], max_loop[可选], output_structure[可选])`
- **输入：** Tm 必填 2–36 nt canonical DNA；heterodimer 必填两条 1–60 nt canonical DNA；可选 Na⁺/K⁺、Mg²⁺、dNTP、链浓度和结构选项。
- **示例代码：**

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

- **当前本地示例结果：**

```text
结果随用户安装的 Primer3 版本、参数文件和条件变化
```

- **限制：** DNAKit 不会自动安装、下载或从 `PATH` 搜索 Primer3；数值依赖实际版本、参数文件和显式条件。`ntthal` CLI 不暴露本 adapter 可安全映射的 DMSO/甲酰胺选项，因此结构计算在这些值非零时会拒绝。adapter 不接受用户预设 alignment 或任意修饰碱基参数。

!!! warning "当前不支持任意修饰碱基热力学"
    修饰碱基需要与具体化学结构和相邻环境匹配的参数集；当前没有公开调用入口，不会用 canonical 参数静默替代。
