# 理化性质

集中说明 DNA 序列的理化性质，以及双链热力学、碱基堆积、稳定性、发卡和二聚体形成倾向。

## 1) THERMO-001 分子量

- **作用：** 根据碱基组成、单双链类型和末端设置估算 DNA 的理论分子量，返回 Da 数值，用于摩尔浓度、质量浓度及样品用量换算。
- **API：** `dnakit.thermodynamics.molecular_weight(sequence[必须], strand[可选], five_prime_phosphorylated[可选], max_sequence_length[可选])`
- **输入：** 必填线性 canonical `DNASequence`；可选单/双链和 5′ 磷酸化状态。
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

- **限制：** 采用明确端基约定的理论近似；不支持修饰碱基，双链按完全互补链计算。

## 2) THERMO-014 260 nm 消光系数

- **作用：** 根据单链 DNA 序列计算 260 nm 理论摩尔消光系数 `ε260`，供 Beer–Lambert 定律换算核酸浓度和吸光度。
- **API：** `dnakit.thermodynamics.extinction_coefficient_260nm(sequence[必须], max_sequence_length[可选])`
- **输入：** 必填 1–1,000,000 nt 的线性、无 Gap、A/C/G/T、单链且未修饰的 `DNASequence`。
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

- **计算口径：** 对长度至少为 2 的序列，使用“所有相邻二核苷酸系数之和，减去所有内部单碱基系数之和”；单核苷酸直接使用其单碱基系数。参数表参考条件为 25 °C、pH 7，单位 `M⁻¹·cm⁻¹`，等价于 `L·mol⁻¹·cm⁻¹`。`ACGT` 的 40,300 结果与公开算例一致。
- **限制：** 这是由序列计算的理论 ε260，不是仪器测得的 A260。A260、光程和该系数可按 Beer–Lambert 定律进一步计算浓度。双链杂交产生的 hypochromicity、荧光染料及其他化学修饰不在当前模型内，不能直接套用本结果。
- **依据：** [IDT 的寡核苷酸定量说明](https://sg.idtdna.com/page/support-and-education/decoded-plus/oligo-quantification-getting-it-right)、[ACGT 公开算例与参数表](https://www.sigmaaldrich.com/US/en/technical-documents/technical-article/genomics/pcr/quantitation-of-oligos)、Warshaw–Tinoco 1966（DOI [`10.1016/0022-2836(66)90115-X`](https://doi.org/10.1016/0022-2836(66)90115-X)）和 Cantor–Warshaw–Shapiro 1970（DOI [`10.1002/bip.1970.360090909`](https://doi.org/10.1002/bip.1970.360090909)）。

## 3) THERMO-002 熔解温度Tm

- **作用：** 根据序列长度、组成或最近邻模型及实验条件计算 DNA 熔解温度 Tm，用于比较双链稳定性和设计寡核苷酸实验条件。
- **API：** `dnakit.thermodynamics.melting_temperature(sequence[必须], method[可选], conditions[可选], config[可选])`
- **输入：** 必填 canonical `DNASequence`；可选 `method`、`ThermodynamicConditions` 和 NN 配置。
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

- **限制：** Wallace 仅接受 2–13 nt 且不建模溶液条件；内部 NN 仅接受 2–60 nt、线性、无 Gap、A/C/G/T 和 Na⁺+K⁺ 总单价盐。`melting_temperature()` 不接收 adapter；Mg²+、dNTP 等条件需用户单独安装 Primer3，并显式构造带 `oligotm_path` 的 `Primer3CLIAdapter.tm()`。

## 4) THERMO-003 盐浓度修正

- **作用：** 根据 Na⁺、K⁺ 等单价离子浓度对 Tm 进行经验修正，返回指定盐环境下的温度变化，便于比较不同缓冲条件。
- **API：** `dnakit.thermodynamics.salt_correction(sequence_length[必须], conditions[可选])`
- **输入：** 必填序列长度；可选包含 Na⁺、K⁺ 浓度的 `ThermodynamicConditions`。
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

- **限制：** 内部模型只修正 Na⁺+K⁺ 总单价盐；Mg²+ 和 dNTP 不会被静默折算或混入。

## 5) THERMO-012 局部熔解特征

- **作用：** 沿序列滑动窗口计算局部 Tm，返回每个窗口的位置、温度及整体范围，用于发现稳定性异常的局部区域。
- **API：** `dnakit.thermodynamics.window_tm(sequence[必须], window_size[必须], step[可选], method[可选], conditions[可选], config[可选], max_windows[可选])`
- **输入：** 必填序列和窗口大小；可选步长、方法、条件、NN 配置和窗口上限。
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

- **限制：** 这是窗口 Tm 数列，不是实验 melting profile；每个窗口仍受所选 Wallace/NN 适用域限制。

## 6) `EVAL-013` 合成风险

- **作用：** 按 GC 极端、同碱基连续、重复和已知规则逐项筛查序列，返回命中位置、风险项和透明分数，用于合成前的规则型预检查。
- **API：** `dnakit.evaluation.evaluate_synthesis_risk(value[必须], config[可选])`；`config` 使用 `dnakit.evaluation.SynthesisRiskConfig`。
- **输入：** 线性、无显式 Gap 的 DNA；可选阈值与资源上限。
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

- **限制：** 不是供应商接单规则、实验合成成功率或真实折叠预测。

## 热力学性质

计算 DNA 双链的热力学参数、碱基堆积、稳定性、发卡结构和二聚体形成倾向。

内部模型和可选外部后端的适用域不同；所有条件浓度单位均为 mol/L。

### 1) THERMO-004 热力学参数

- **作用：** 使用最近邻参数和显式实验条件计算完全互补双链形成时的 ΔH、ΔS、ΔG 和 Tm，用于定量比较候选双链的热力学稳定性。
- **API：** `dnakit.thermodynamics.nearest_neighbor(sequence[必须], complement[可选], conditions[可选], config[可选])`
- **输入：** 必填一条 canonical 序列；可选完全互补链、条件和 NN 参数集。
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

- **限制：** 仅支持完全 Watson–Crick 互补的内部参数域；mismatch、dangling end 和修饰不支持。

### 2) THERMO-005 Nearest-neighbor

- **作用：** 把双链拆成相邻碱基对步骤，分别计算每一步对 ΔH、ΔS 和 ΔG 的贡献，并汇总整条序列，便于定位稳定或不稳定片段。
- **API：** `dnakit.thermodynamics.nearest_neighbor(sequence[必须], complement[可选], conditions[可选], config[可选])`
- **输入：** 必填 2–60 nt canonical 线性序列；可选完全互补链、条件和 `NearestNeighborConfig`。
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

- **限制：** 参数集当前固定为 `santalucia1998-v1`；不是 Primer3、MELTING 或 NUPACK 的通用替代。

### 3) THERMO-006 Duplex stability

- **作用：** 检查两条 DNA 是否满足完整 Watson–Crick 互补关系；通过检查后计算双链热力学结果，用于验证配对序列及其稳定性。
- **API：** `dnakit.thermodynamics.duplex_stability(sequence_a[必须], sequence_b[必须], conditions[可选], config[可选], backend[可选], adapter[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填序列 A、B；可选条件、`native`/`primer3-cli` 后端、adapter、max loop 和结构输出。
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

- **限制：** 默认 `native` 只接受完全反向互补。canonical mismatch/dangling heterodimer 需用户单独安装 Primer3、构造含显式 `ntthal_path` 的 `Primer3CLIAdapter`，并选择 `backend="primer3-cli"`。接口不会从 `PATH` 自行构造可用后端。不支持修饰或用户预设 alignment。

### 4) THERMO-007 碱基堆积

- **作用：** 单独列出每个相邻碱基堆积步骤对焓、熵和自由能的贡献，用于解释整条双链热力学结果来自哪些局部序列。
- **API：** `dnakit.thermodynamics.stacking_interactions(sequence[必须], temperature_celsius[可选], config[可选])`
- **输入：** 必填 canonical `DNASequence`；可选温度和 NN 配置。
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

- **限制：** 只报告 stacking 项，不含 initiation、terminal、symmetry 或 salt 修正，不能单独解释为完整 duplex ΔG。

### 5) THERMO-008 Hairpin

- **作用：** 调用显式结构后端预测单条序列是否形成发卡，并返回结构是否存在、Tm 和 ΔG 等指标，用于筛查可能影响引物或寡核苷酸使用的自折叠。
- **API：** `dnakit.thermodynamics.probe_primer3(ntthal_path[可选], thermodynamic_parameters_path[可选])`、`dnakit.thermodynamics.Primer3CLIAdapter.hairpin(sequence[必须], conditions[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填 1–60 nt canonical 序列；可选条件、`max_loop=1..30` 和结构输出。
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

- **限制：** DNAKit 不自动安装、下载或从 `PATH` 搜索 Primer3；许可提示为 `GPL-2.0-or-later`，应以实际安装版本为准。当前自动测试验证 CLI 协议和解析边界，不是 native、NUPACK 或真实 Primer3 科学差分。

### 6) THERMO-009 Self-dimer

- **作用：** 评估同一条序列的两个分子是否可能形成自二聚体，返回预测结构、Tm 和 ΔG，用于筛查引物自身配对风险。
- **API：** `dnakit.thermodynamics.Primer3CLIAdapter.self_dimer(sequence[必须], conditions[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填 1–60 nt canonical 序列；可选条件、max loop 和结构输出。
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

- **限制：** 仅在调用方合法安装并显式调用 `ntthal` 时执行；不会自动选择后端，也不是实验二聚体测量。

### 7) THERMO-010 Heterodimer

- **作用：** 评估两条不同序列之间是否可能形成异二聚体，返回预测结构、Tm 和 ΔG，用于检查引物对或寡核苷酸之间的非目标配对。
- **API：** `dnakit.thermodynamics.Primer3CLIAdapter.heterodimer(sequence_a[必须], sequence_b[必须], conditions[可选], max_loop[可选], output_structure[可选])`
- **输入：** 必填两条 1–60 nt canonical 序列；可选条件、max loop 和结构输出。
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

- **限制：** 仅支持 Primer3 adapter 接受的 canonical 输入和参数域；结果不是 NUPACK 多复合体平衡或实验稳定性。
