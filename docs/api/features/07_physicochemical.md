# 理化性质

集中说明 DNA 序列的理化性质，以及双链热力学、碱基堆积、稳定性、发卡和二聚体形成倾向。

各项计算的论文来源、内部规则和适用边界见 [FAQ：理化性质的计算依据和参考文献](../../faq.md#physicochemical-references)。

内部模型和可选外部后端的适用域不同；所有条件浓度单位均为 mol/L。

## 1) 分子量 {#1-thermo-001}

- **作用：** 根据碱基组成、单双链类型和末端设置估算 DNA 的理论分子量，返回 Da 数值，用于摩尔浓度、质量浓度及样品用量换算。
- **计算方法：** 对每条链的 A/C/G/T 无水脱氧核苷酸残基质量求和，减去未磷酸化末端修正 `61.96 Da`；若 5′ 磷酸化则增加 `79.0 Da`，双链模式再加上完整反向互补链的质量。
- **API：** `dnakit.thermodynamics.molecular_weight(sequence[必须], strand[可选], five_prime_phosphorylated[可选], max_sequence_length[可选])`
- **输入：** 必填线性 canonical `DNASequence`；可选单/双链和 5′ 磷酸化状态。
- **输出：** 返回 `MolecularWeightResult`，示例展示 `value_dalton`、`value_kilodalton`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import molecular_weight

result = molecular_weight(DNASequence("ACGT"))
print(result.value_dalton, result.value_kilodalton)
```

- **示例结果：**

```text
1173.84 1.17384
```

## 2) 260 nm 消光系数 {#2-thermo-014-260-nm}

- **作用：** 根据单链 DNA 序列计算 260 nm 理论摩尔消光系数 `ε260`，供 Beer–Lambert 定律换算核酸浓度和吸光度。
- **计算方法：** 长度至少为 2 时，将所有相邻二核苷酸的公开消光系数求和，再减去所有内部单碱基系数；单核苷酸直接使用对应的单碱基系数。
- **API：** `dnakit.thermodynamics.extinction_coefficient_260nm(sequence[必须], max_sequence_length[可选])`
- **输入：** 必填 1–1,000,000 nt 的线性、无 Gap、A/C/G/T、单链且未修饰的 `DNASequence`。
- **输出：** 返回 `ExtinctionCoefficientResult`，示例展示 `value_m_inverse_cm_inverse`、`wavelength_nm`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import extinction_coefficient_260nm

result = extinction_coefficient_260nm(DNASequence("ACGT"))
print(result.value_m_inverse_cm_inverse, result.wavelength_nm)
```

- **示例结果：**

```text
40300.0 260
```


## 3) 熔解温度Tm {#3-thermo-002-tm}

- **作用：** 根据序列长度、组成或最近邻模型及实验条件计算 DNA 熔解温度 Tm，用于比较双链稳定性和设计寡核苷酸实验条件。
- **计算方法：** `wallace` 方法使用 `2 × (A+T) + 4 × (G+C)`；`nearest_neighbor` 方法使用 SantaLucia 1998 相邻碱基堆积、末端、对称性、盐浓度和链浓度参数计算 Tm。
- **API：** `dnakit.thermodynamics.melting_temperature(sequence[必须], method[可选], conditions[可选], config[可选])`
- **输入：** 必填 canonical `DNASequence`；可选 `method`、`ThermodynamicConditions` 和 NN 配置。
- **输出：** 返回 `MeltingTemperatureResult`，示例展示 `tm_celsius`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import melting_temperature

result = melting_temperature(DNASequence("AACG"), method="wallace")
print(result.tm_celsius)
```

- **示例结果：**

```text
12.0
```

## 4) 盐浓度修正 {#4-thermo-003}

- **作用：** 根据 Na⁺、K⁺ 等单价离子浓度对 Tm 进行经验修正，返回指定盐环境下的温度变化，便于比较不同缓冲条件。
- **计算方法：** 使用 SantaLucia 1998 单价盐熵修正式 `ΔS_salt = 0.368 × (N - 1) × ln([Na⁺] + [K⁺])`，其中 `N` 为序列长度、浓度单位为 mol/L。
- **API：** `dnakit.thermodynamics.salt_correction(sequence_length[必须], conditions[可选])`
- **输入：** 必填序列长度；可选包含 Na⁺、K⁺ 浓度的 `ThermodynamicConditions`。
- **输出：** 返回 `SaltCorrectionResult`，示例展示 `delta_s_cal_per_k_mol`。
- **示例代码：**

```python
from dnakit.thermodynamics import ThermodynamicConditions, salt_correction

conditions = ThermodynamicConditions(sodium_molar=0.05)
result = salt_correction(10, conditions=conditions)
print(round(result.delta_s_cal_per_k_mol, 3))
```

- **示例结果：**

```text
-9.922
```

## 5) 局部熔解特征 {#5-thermo-012}

- **作用：** 沿序列滑动窗口计算局部 Tm，返回每个窗口的位置、温度及整体范围，用于发现稳定性异常的局部区域。
- **计算方法：** 按 `window_size` 和 `step` 生成固定窗口，对每个窗口重复调用所选的 Wallace 或 SantaLucia nearest-neighbor Tm 方法，并汇总局部最小值、最大值和位置。
- **API：** `dnakit.thermodynamics.window_tm(sequence[必须], window_size[必须], step[可选], method[可选], conditions[可选], config[可选], max_windows[可选])`
- **输入：** 必填序列和窗口大小；可选步长、方法、条件、NN 配置和窗口上限。
- **输出：** 返回 `WindowTmResult`，示例展示 `windows`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import window_tm

result = window_tm(DNASequence("AACGTT"), 4, step=2, method="wallace")
print([(item.start, item.end, item.tm_celsius) for item in result.windows])
```

- **示例结果：**

```text
[(0, 4, 12.0), (2, 6, 12.0)]
```

## 6) 合成风险 {#6-eval-013}

- **作用：** 按 GC 极端、同碱基连续、重复和已知规则逐项筛查序列，返回命中位置、风险项和透明分数，用于合成前的规则型预检查。
- **计算方法：** 分别计算全局 GC、局部 GC 异常窗口、最长同碱基连续、串联重复和倒置重复五个风险分量，截断到 `[0, 1]` 后等权求平均；分数 `<0.2`、`0.2–<0.5`、`≥0.5` 依次标为低、中、高风险。
- **API：** `dnakit.evaluation.evaluate_synthesis_risk(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.SynthesisRiskConfig`。
- **输入：** 线性、无显式 Gap 的 DNA；可选阈值与资源上限。
- **输出：** 返回 `EvaluationReport`，整体指标在 `metrics`，逐条结果在 `entries`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.evaluation import evaluate_synthesis_risk

report = evaluate_synthesis_risk(DNASequence("G" * 40 + "AT" * 10))
entry = report.entries[0]
print(entry.metrics["risk_level"])
print(entry.metrics["risk_score"])
```

- **示例结果：**

```text
medium
0.38666666666666666
```

<span id="_2"></span>

## 7) 热力学参数 {#1-thermo-004}

- **作用：** 使用最近邻参数和显式实验条件计算完全互补双链形成时的 ΔH、ΔS、ΔG 和 Tm，用于定量比较候选双链的热力学稳定性。
- **计算方法：** 将 SantaLucia 1998 的相邻堆积、末端起始、对称性和单价盐贡献分别求和得到 ΔH 与 ΔS，再用 `ΔG = ΔH - TΔS` 及浓度修正的热力学方程计算 Tm。
- **API：** `dnakit.thermodynamics.nearest_neighbor(sequence[必须], complement[可选], conditions[可选], config[可选])`
- **输入：** 必填一条 canonical 序列；可选完全互补链、条件和 NN 参数集。
- **输出：** 返回 `NearestNeighborResult`，示例展示 `delta_h_kcal_per_mol`、`delta_s_cal_per_k_mol`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, nearest_neighbor

conditions = ThermodynamicConditions(sodium_molar=1.0, strand_concentration_molar=1e-6)
result = nearest_neighbor(DNASequence("GTGCAT"), conditions=conditions)
print(result.delta_h_kcal_per_mol, result.delta_s_cal_per_k_mol)
```

- **示例结果：**

```text
-40.0 -111.3
```

## 8) Nearest-neighbor {#2-thermo-005-nearest-neighbor}

- **作用：** 把双链拆成相邻碱基对步骤，分别计算每一步对 ΔH、ΔS 和 ΔG 的贡献，并汇总整条序列，便于定位稳定或不稳定片段。
- **计算方法：** 将序列拆成所有相邻二核苷酸步骤，逐项查询 SantaLucia 1998 ΔH/ΔS 参数表，再加入末端、对称性和盐修正；该项与热力学参数使用同一个 `nearest_neighbor()` 结果，只是强调逐步骤明细。
- **API：** `dnakit.thermodynamics.nearest_neighbor(sequence[必须], complement[可选], conditions[可选], config[可选])`
- **输入：** 必填 2–60 nt canonical 线性序列；可选完全互补链、条件和 `NearestNeighborConfig`。
- **输出：** 返回 `NearestNeighborResult`，示例展示 `stacking_steps`、`tm_celsius`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import ThermodynamicConditions, nearest_neighbor

result = nearest_neighbor(
    DNASequence("GTGCAT"),
    conditions=ThermodynamicConditions(sodium_molar=1.0, strand_concentration_molar=1e-6),
)
print([step.top_5to3 for step in result.stacking_steps], round(result.tm_celsius, 3))
```

- **示例结果：**

```text
['GT', 'TG', 'GC', 'CA', 'AT'] 9.517
```

## 9) Duplex stability {#3-thermo-006-duplex-stability}

- **作用：** 检查两条 DNA 是否满足完整 Watson–Crick 互补关系；通过检查后计算双链热力学结果，用于验证配对序列及其稳定性。
- **计算方法：** `native` 后端要求两条链完全反向互补并使用 SantaLucia nearest-neighbor 结果，以 `Tm > 设置温度` 判为稳定；显式选择 `primer3-cli` 时调用 Primer3 `ntthal` 估算允许 mismatch 和 dangling end 的异二聚体结构。
- **API：** `dnakit.thermodynamics.duplex_stability(sequence_a[必须], sequence_b[必须], conditions[可选], config[可选], backend[可选], adapter[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填序列 A、B；可选条件、`native`/`primer3-cli` 后端、adapter、max loop 和结构输出。
- **输出：** 返回 `DuplexStabilityResult`，示例展示 `fully_complementary`、`stable_at_temperature`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import duplex_stability

result = duplex_stability(DNASequence("GTGCAT"), DNASequence("ATGCAC"))
print(result.fully_complementary, result.stable_at_temperature)
```

- **示例结果：**

```text
True False
```

## 10) 碱基堆积 {#4-thermo-007}

- **作用：** 单独列出每个相邻碱基堆积步骤对焓、熵和自由能的贡献，用于解释整条双链热力学结果来自哪些局部序列。
- **计算方法：** 对每个相邻二核苷酸查询 SantaLucia 1998 的 ΔH 和 ΔS 参数，并按指定温度计算 `ΔG = ΔH - TΔS` 后求和；该接口不加入末端、对称性和盐修正。
- **API：** `dnakit.thermodynamics.stacking_interactions(sequence[必须], temperature_celsius[可选], config[可选])`
- **输入：** 必填 canonical `DNASequence`；可选温度和 NN 配置。
- **输出：** 返回 `StackingResult`，示例展示 `steps`、`total_delta_h_kcal_per_mol`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.thermodynamics import stacking_interactions

result = stacking_interactions(DNASequence("AA"))
print(result.steps[0].top_5to3, result.total_delta_h_kcal_per_mol)
```

- **示例结果：**

```text
AA -7.9
```

## 11) Hairpin {#5-thermo-008-hairpin}

- **作用：** 调用显式结构后端预测单条序列是否形成发卡，并返回结构是否存在、Tm 和 ΔG 等指标，用于筛查可能影响引物或寡核苷酸使用的自折叠。
- **计算方法：** 调用用户独立安装的 Primer3 `ntthal`，以 `HAIRPIN` 模式在给定盐浓度、链浓度、温度和 `max_loop` 条件下搜索单分子发卡结构并解析其 Tm、ΔG、ΔH 和 ΔS。
- **API：** `dnakit.thermodynamics.probe_primer3(ntthal_path[可选], thermodynamic_parameters_path[可选])`、`dnakit.thermodynamics.Primer3CLIAdapter.hairpin(sequence[必须], conditions[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填 1–60 nt canonical 序列；可选条件、`max_loop=1..30` 和结构输出。
- **输出：** `probe_primer3` 返回后端可用性信息；发卡分析返回结构是否形成及 Tm 等热力学结果。
- **示例代码：**

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

- **示例结果：**

```text
结果随用户安装的 Primer3 版本、参数文件和条件变化
```

## 12) Self-dimer {#6-thermo-009-self-dimer}

- **作用：** 评估同一条序列的两个分子是否可能形成自二聚体，返回预测结构、Tm 和 ΔG，用于筛查引物自身配对风险。
- **计算方法：** 将同一条序列作为两条输入调用用户安装的 Primer3 `ntthal` `ANY` 模式，在指定实验条件下搜索自二聚体热力学结构并解析相关参数。
- **API：** `dnakit.thermodynamics.Primer3CLIAdapter.self_dimer(sequence[必须], conditions[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填 1–60 nt canonical 序列；可选条件、max loop 和结构输出。
- **输出：** 返回 `Primer3ThermodynamicResult`，示例展示 `structure_found`、`delta_g_kcal_per_mol`。
- **示例代码：**

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

- **示例结果：**

```text
结果随用户安装的 Primer3 版本、参数文件和条件变化
```

## 13) Heterodimer {#7-thermo-010-heterodimer}

- **作用：** 评估两条不同序列之间是否可能形成异二聚体，返回预测结构、Tm 和 ΔG，用于检查引物对或寡核苷酸之间的非目标配对。
- **计算方法：** 将两条不同序列输入用户安装的 Primer3 `ntthal` `ANY` 模式，在指定实验条件下搜索异二聚体热力学结构并解析相关参数。
- **API：** `dnakit.thermodynamics.Primer3CLIAdapter.heterodimer(sequence_a[必须], sequence_b[必须], conditions[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填两条 1–60 nt canonical 序列；可选条件、max loop 和结构输出。
- **输出：** 返回 `Primer3ThermodynamicResult`，示例展示 `structure_found`、`delta_g_kcal_per_mol`。
- **示例代码：**

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

- **示例结果：**

```text
结果随用户安装的 Primer3 版本、参数文件和条件变化
```
