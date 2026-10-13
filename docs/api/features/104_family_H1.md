# H1 家族变体复核：dadi、tskit 与 scikit-allel

本页记录第 H1 组对台账中 `family-variant` 行的复核。这些行原先被判定为已包装函数的别名、再导出或由参数选择的分支；这里逐行核对该说法是否成立，并把确实有独立 DNA 级含义与输出的函数加以包装。新增函数位于 `src/dnakit/population/dadi_family_functions.py`、`src/dnakit/population/tskit_family_functions.py` 与 `src/dnakit/variants/allel_family_functions.py`。它们都已从 `dnakit.population` 与 `dnakit.variants` 导出，示例从这些包导入。

- 新增 26 个函数：dadi 6 个、tskit 2 个、scikit-allel 18 个，每个都在下面有带真实输出的示例。
- 测试位于 `others/tests/unit/test_family_H1.py`。每个新函数与上游调用或手工计算对照；保留为 `family-variant` 的行中有一批按声称的参数调用已有包装函数，并与上游变体的结果对照。
- 示例中的输出由实际运行得到，文档检查器应报告 0 个不一致。

## 复核结果

| 包 | 复核前 family-variant | duplicate-of | family-variant（已核实） | wrapped | 其他状态 |
| --- | --- | --- | --- | --- | --- |
| dadi | 121 | 1 | 108 | 10 | broken-upstream 2 |
| tskit | 11 | 0 | 9 | 2 | - |
| allel | 28 | 1 | 0 | 25 | broken-upstream 1; needs-binary 1 |

合计 160 行：保留 117 行 `family-variant`，37 行改为 `wrapped`，2 行 `duplicate-of`，其余 4 行改为 `broken-upstream` 或 `needs-binary`。

### 保留为 family-variant 的行（已核实）

保留的前提是：已有包装函数的签名里确有该参数（用 `inspect.signature` 核对），且包装函数内部调用的正是该行对应的上游对象。

| 包装函数与参数 | 行数 | 对应的上游函数 |
| --- | --- | --- |
| `expected_sfs_1d(model=...)` | 5 | `Demographics1D.growth`、`bottlegrowth`、`bottlegrowth_1d`、`three_epoch`、`three_epoch_inbreeding` |
| `expected_sfs_2d(model=...)` | 42 | `Demographics2D` 的其余双种群模型（`no_mig`、`IM`、`founder_asym` 等） |
| `expected_sfs_3d(model=...)` | 33 | `Demographics3D` 的其余三种群模型（`out_of_africa` 等） |
| `expected_sfs_selection(model=...)` | 20 | `DFE.DemogSelModels` 的其余选择模型 |
| `ms_core_command_2d(model='IM' / 'IM_pre')` | 2 | `Demographics2D.IM_mscore`、`IM_pre_mscore` |
| `nested_model_statistic(test='wald' / 'score')` | 2 | `Godambe.Wald_stat`、`score_stat` |
| `fit_demographic_model(method='nelder-mead' / 'powell')` | 2 | `Inference.optimize_log_fmin`、`optimize_log_powell` |
| `sfs_residuals(residual='poisson')` | 1 | `Inference.linear_Poisson_residual` |
| `read_vcf_data(parser='cyvcf2')` | 1 | `Misc.make_data_dict_vcf_cyvcf2` |
| `read_text_table(table=...)` | 6 | `tskit.parse_individuals`、`parse_edges`、`parse_sites`、`parse_mutations`、`parse_populations`、`parse_migrations` |
| `random_tree_newick(shape='star' / 'comb' / 'balanced')` | 3 | `tskit.combinatorics.generate_star`、`generate_comb`、`generate_balanced` |

dadi 的 1D、2D、3D 与选择模型共 100 行都用实际调用核实：每个模型名都在包装函数的固定列表里，`getattr` 取到的对象就是该行的对象，并且用小样本量和小网格实际运行成功。

### 被证伪的说法

- dadi 的 5 个似然函数（`ll_multinom`、`minus_ll`、`minus_ll_multinom`、`ll_per_bin`、`ll_multinom_per_bin`）与 `optimally_scaled_sfs` 被说成是 `calculate_sfs_log_likelihood`、`estimate_sfs_theta` 的变体，台账原文也写着 “not exposed”。这两个包装函数只对中性模型和一维谱起作用，没有任何参数能到达它们，因此改为 `wrapped`（`calculate_model_log_likelihood`、`scale_model_sfs_to_data`）。
- dadi 的 `optimize`、`optimize_cons`、`optimize_grid` 与 `opt` 被说成 `fit_demographic_model` 的变体。该函数的 `method` 只有四个对数尺度的取值，没有约束、网格或 NLopt 选项，因此这四行改为 `wrapped`。`NLopt_mod.opt` 与 `Inference.opt` 是同一个对象（`is` 为真），改为 `duplicate-of`。
- `optimize_log_lbfgsb` 被说成可由 `fit_demographic_model(method='lbfgsb')` 到达，但 dadi 2.4.4 在调用 `scipy.optimize.fmin_l_bfgs_b` 时传入 `iprint`，SciPy 1.18.1 已不接受，上游直接抛出 TypeError（`optimize_lbfgsb` 同样）。两行都改为 `broken-upstream`。需要注意的是 `fit_demographic_model` 的文档仍把 `lbfgsb` 列为可用方法，在当前环境里它总是得到 ConfigurationError。
- scikit-allel 的 28 行中，只有 `average_patterson_d` 的说法成立，但它不是参数变体，而是 `blockwise_patterson_d` 的同一个对象（`duplicate-of`）。其余各行的窗口形式、逐变异形式、联合谱、折叠与缩放形式、刀切误差形式、随机化 PCA、分步定相和字段表格式，都没有任何已有包装函数的参数能到达，其中 25 行改为 `wrapped`。例如 `calculate_fst` 只有 `method` 一个参数，只返回平均值。
- `allel.roh_poissonhmm` 需要可选包 protopunica，完整环境里没有，也无法在此安装，改为 `needs-binary`。`allel.write_vcf` 上游自述为初步实现，实测写出的 VCF 不合规（ALT 为 `G,`，缺失 QUAL 写成 `nan`，FILTER 全部为 PASS），改为 `broken-upstream`。
- tskit 的 `write_nexus` 被说成可由 `export_sample_sequences(file_format='nexus')` 到达，但该包装函数固定 `include_trees=False`，不能写出 TREES 块；`tree_count_topologies` 被说成 `count_sample_topologies` 的单树版本，但后者没有 `tree_index` 参数，会计算每一棵树。两行都改为 `wrapped`。

### 复核中发现的其他问题（不在本组行内）

- `allel.fold_sfs` 对偶数个染色体给出错误结果：对 `[1, 2, 7, 14, 11, 5, 0]`（n=6）返回 `[6, 13, 21]`，正确的折叠谱是 `[1, 7, 18, 14]`（类别 i 与 n-i 相加，`sfs_folded` 给出这个结果，dadi 的折叠对未掩盖的类别 1 到 3 也给出 7、18、14）。台账中该行被标为与 dadi 折叠 `duplicate-of`，对偶数 n 并不成立。
- `expected_sfs_selection` 等现有包装函数只把 ValueError、TypeError、IndexError、KeyError 转为 ConfigurationError。dadi 的 `IM_pre_sel` 在某些参数下抛出 ZeroDivisionError，会直接传给调用者。

## 通用约定

- 每个公开函数的 `verbose` 都是仅关键字参数，默认 `True`。
- 返回值是可 JSON 序列化的字典，包含 `status`、`operation`、结果字段与 `provenance`。
- 统计量无法计算时输出 `null`；dadi 掩盖的条目也输出 `null`。
- dadi、tskit、scikit-allel 只在函数运行时导入；缺少依赖时抛出 `BackendUnavailableError`。
- 输入不合法时抛出 `ConfigurationError`。输出文件不会覆盖已存在的文件。
- 随机化 PCA 的随机数只通过 `seed` 固定，不修改任何全局随机状态。

## 1) `calculate_model_log_likelihood`

- **作用：** 给定任意模型谱，计算观测计数的对数似然。
- **处理方法：** 把 `model` 与 `data` 转为 dadi 谱（角点被掩盖），再调用 `Inference.ll`；`multinom=True` 用 `ll_multinom`（先把模型缩放到最优 theta）；`negate=True` 用 `minus_ll` 或 `minus_ll_multinom`；`per_bin=True` 用 `ll_per_bin` 或 `ll_multinom_per_bin`，返回每个条目的值。`negate` 与 `per_bin` 不能同时使用。已有的 `calculate_sfs_log_likelihood` 只给中性模型的 Poisson 总和，没有这些选项。
- **API：** `calculate_model_log_likelihood(model[必须], data[必须], *, multinom[可选], per_bin[可选], negate[可选], verbose[可选])`
- **输入：** `model` 与 `data` 是形状相同的嵌套列表（一到三个种群，每个维度至少 3 项）；`data` 必须非负。
- **输出：** 字典，含 `dadi_function`（实际调用的 dadi 函数）、`multinom`、`per_bin`、`negate`、`shape`，以及 `log_likelihood`（标量）或 `log_likelihood_per_bin`（嵌套列表，掩盖项为 null）。

**示例代码：**

```python
from dnakit.population import calculate_model_log_likelihood

model = [0.0, 10.0, 5.0, 3.3, 2.5, 2.0, 1.7, 1.4, 1.2, 1.0, 0.0]
data = [0.0, 12.0, 4.0, 4.0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 0.0]
plain = calculate_model_log_likelihood(model, data, verbose=False)
multi = calculate_model_log_likelihood(model, data, multinom=True, verbose=False)
negated = calculate_model_log_likelihood(model, data, negate=True, verbose=False)
per_bin = calculate_model_log_likelihood(model, data, per_bin=True, verbose=False)
print(plain["dadi_function"], round(plain["log_likelihood"], 4))
print(multi["dadi_function"], round(multi["log_likelihood"], 4))
print(negated["dadi_function"], round(negated["log_likelihood"], 4))
print([None if v is None else round(v, 3) for v in per_bin["log_likelihood_per_bin"]])
```

**示例输出：**

```text
ll -13.0735
ll_multinom -13.0593
minus_ll 13.0735
[None, -2.356, -1.74, -1.702, -1.361, -1.307, -1.169, -1.42, -1.018, -1.0, None]
```

## 2) `scale_model_sfs_to_data`

- **作用：** 把模型谱乘以数据的最大似然 theta。
- **处理方法：** 调用 `Inference.optimally_scaled_sfs` 得到缩放后的谱，并用 `Inference.optimal_sfs_scaling` 给出缩放因子 theta（数据总和除以模型总和，只计未被掩盖的条目）。已有的 `estimate_sfs_theta` 只对中性模型给出 theta，不返回缩放后的谱。
- **API：** `scale_model_sfs_to_data(model[必须], data[必须], *, verbose[可选])`
- **输入：** `model` 与 `data` 同上，形状必须相同。
- **输出：** 字典，含 `shape`、`theta` 与 `scaled_spectrum`（嵌套列表，掩盖项为 null）。

**示例代码：**

```python
from dnakit.population import scale_model_sfs_to_data

model = [0.0, 10.0, 5.0, 3.3, 2.5, 2.0, 1.7, 1.4, 1.2, 1.0, 0.0]
data = [0.0, 12.0, 4.0, 4.0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 0.0]
result = scale_model_sfs_to_data(model, data, verbose=False)
print(round(result["theta"], 6))
print([None if v is None else round(v, 3) for v in result["scaled_spectrum"]])
```

**示例输出：**

```text
1.032028
[None, 10.32, 5.16, 3.406, 2.58, 2.064, 1.754, 1.445, 1.238, 1.032, None]
```

## 3) `fit_demographic_model_natural_scale`

- **作用：** 用 BFGS 在原始（非对数）参数尺度上拟合人口历史模型。
- **处理方法：** 调用 `Inference.optimize`。已有的 `fit_demographic_model` 的所有 `method` 取值都在对数参数上优化，没有参数能到达这个函数。拟合后用 `ll_multinom`（或 `ll`）给出对数似然。
- **API：** `fit_demographic_model_natural_scale(data[必须], *, model[必须], p0[必须], lower_bound[可选], upper_bound[可选], multinom[可选], grid_points[可选], maxiter[可选], verbose[可选])`
- **输入：** `data` 是一或两个种群的观测谱（嵌套列表）；`model` 是该维度模型族中的模型名；`p0` 是正的起始参数；`lower_bound`、`upper_bound` 与 `p0` 等长。`multinom=True` 时 theta 在每一步重新缩放，不属于 `p0`。
- **输出：** 字典，含 `params`、`log_likelihood`、`function_calls`、`converged`（BFGS 的 warnflag 为 0）以及所用选项。

**示例代码：**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_natural_scale

result = fit_demographic_model_natural_scale(
    data,
    model="two_epoch",
    p0=[1.5, 0.7],
    lower_bound=[1e-3, 1e-3],
    upper_bound=[100, 10],
    maxiter=40,
    verbose=False,
)
print([round(v, 2) for v in result["params"]], round(result["log_likelihood"], 2))
```

**示例输出：**

```text
[2.0, 0.5] -33.55
```

## 4) `fit_demographic_model_constrained`

- **作用：** 在线性参数约束下拟合人口历史模型。
- **处理方法：** 调用 `Inference.optimize_cons`（SLSQP）。每个约束是 `{"coefficients": [...], "offset": b}`，表示 `c·参数 + b`；等式约束在解处必须为 0，不等式约束必须不小于 0。例如 `nu1 = 2·nu2` 写作 `{"coefficients": [1, -2]}`。已有的 `fit_demographic_model` 没有约束参数。
- **API：** `fit_demographic_model_constrained(data[必须], *, model[必须], p0[必须], lower_bound[可选], upper_bound[可选], equality_constraints[可选], inequality_constraints[可选], multinom[可选], grid_points[可选], maxiter[可选], verbose[可选])`
- **输入：** `data` 是一或两个种群的观测谱（嵌套列表）；`model` 是该维度模型族中的模型名；`p0` 是正的起始参数；`lower_bound`、`upper_bound` 与 `p0` 等长。`multinom=True` 时 theta 在每一步重新缩放，不属于 `p0`。 `equality_constraints` 与 `inequality_constraints` 是约束映射的列表，至少给出一个。
- **输出：** 字典，含 `params`、`log_likelihood`、`equality_values` 与 `inequality_values`（解处各约束的取值）、`iterations`、`converged`。

**示例代码：**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_constrained

result = fit_demographic_model_constrained(
    data,
    model="two_epoch",
    p0=[1.5, 0.7],
    lower_bound=[1e-3, 1e-3],
    upper_bound=[100, 10],
    equality_constraints=[{"coefficients": [1, -4]}],
    inequality_constraints=[{"coefficients": [1, -1]}],
    maxiter=40,
    verbose=False,
)
nu, time = result["params"]
print(round(nu, 2), round(time, 2), round(nu - 4 * time, 6), result["converged"])
```

**示例输出：**

```text
2.0 0.5 0.0 True
```

## 5) `fit_demographic_model_grid`

- **作用：** 在参数网格上穷举搜索，拟合人口历史模型。
- **处理方法：** 调用 `Inference.optimize_grid`。`parameter_ranges` 对每个模型参数给出 `[start, stop, num_points]`，网格取 `start` 到 `stop`（含端点）的等距点，全部组合（至多 2000 个）都会被评估。适合在局部拟合前寻找起始值。已有的 `fit_demographic_model` 没有网格搜索。
- **API：** `fit_demographic_model_grid(data[必须], *, model[必须], parameter_ranges[必须], multinom[可选], grid_points[可选], verbose[可选])`
- **输入：** `data`、`model`、`multinom` 与 `grid_points` 同上；`parameter_ranges` 的条目数必须等于模型参数个数，`start` 必须小于 `stop`，`num_points` 至少 2。
- **输出：** 字典，含 `params`（最佳网格点）、`log_likelihood`、`num_evaluations`、`parameter_ranges` 与所用选项。

**示例代码：**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_grid

result = fit_demographic_model_grid(
    data,
    model="two_epoch",
    parameter_ranges=[[1.0, 3.0, 3], [0.2, 1.0, 5]],
    verbose=False,
)
print([round(v, 2) for v in result["params"]], result["num_evaluations"])
print(round(result["log_likelihood"], 2))
```

**示例输出：**

```text
[2.0, 0.6] 15
-33.77
```

## 6) `fit_demographic_model_nlopt`

- **作用：** 用 NLopt 的无导数算法拟合人口历史模型。
- **处理方法：** 调用 `Inference.opt`（与 `NLopt_mod.opt` 是同一个对象）。`algorithm` 可取 `LN_BOBYQA`（默认）、`LN_COBYLA`、`LN_NELDERMEAD` 或 `LN_SBPLX`；需要 `nlopt` 包。搜索在 `maxeval` 次模型计算后，或对数似然变化小于 `ftol_abs`、参数变化小于 `xtol_abs` 时停止。上游的 `log_opt` 选项会返回错误的参数，因此不开放。
- **API：** `fit_demographic_model_nlopt(data[必须], *, model[必须], p0[必须], lower_bound[可选], upper_bound[可选], algorithm[可选], multinom[可选], grid_points[可选], maxeval[可选], ftol_abs[可选], xtol_abs[可选], verbose[可选])`
- **输入：** `data` 是一或两个种群的观测谱（嵌套列表）；`model` 是该维度模型族中的模型名；`p0` 是正的起始参数；`lower_bound`、`upper_bound` 与 `p0` 等长。`multinom=True` 时 theta 在每一步重新缩放，不属于 `p0`。
- **输出：** 字典，含 `params`、`log_likelihood`、`optimizer_value`（NLopt 报告的最优目标值）、`algorithm` 与所用选项。

**示例代码：**

```python
import dadi
import numpy as np

model = dadi.Demographics1D.two_epoch((2.0, 0.5), (10,), 40)
data = [float(x) for x in np.round(np.ma.filled(model * 1000, 0.0))]
from dnakit.population import fit_demographic_model_nlopt

result = fit_demographic_model_nlopt(
    data,
    model="two_epoch",
    p0=[1.5, 0.7],
    lower_bound=[1e-3, 1e-3],
    upper_bound=[100, 10],
    algorithm="LN_SBPLX",
    maxeval=300,
    verbose=False,
)
print(result["algorithm"], [round(v, 2) for v in result["params"]], round(result["log_likelihood"], 2))
```

**示例输出：**

```text
LN_SBPLX [2.0, 0.5] -33.55
```

## 7) `count_tree_topologies`

- **作用：** 统计一棵指定树中样本集合内嵌的拓扑数。
- **处理方法：** 调用 `tskit.combinatorics.tree_count_topologies`，只对 `tree_index` 指定的那一棵树计算。已有的 `count_sample_topologies` 对每一棵树都计算，没有 `tree_index` 参数；两者对同一棵树给出相同的计数（测试中核对）。
- **API：** `count_tree_topologies(tables[必须], sample_sets[必须], *, tree_index[可选], max_rows[可选], verbose[可选])`
- **输入：** `tables` 是 `build_tree_sequence()` 使用的显式表格映射；`sample_sets` 把名称映射到样本节点 ID 列表，每组至少 2 个样本；`tree_index` 默认 0。
- **输出：** 字典，含 `tree_index`、`interval`（该树的 `[left, right)`）、`sample_sets` 与 `counts`。`counts` 的键是用 `+` 连接的集合名，值是按拓扑秩（`shape`、`label`）排序的 `{shape, label, count}` 列表。

**示例代码：**

```python
TABLES = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 2}, {"time": 3}],
    "edges": [
        {"left": 0, "right": 100, "parent": parent, "child": child}
        for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 20, "ancestral_state": "A"},
        {"position": 30, "ancestral_state": "A"},
    ],
    "mutations": [
        {"site": 0, "node": 0, "derived_state": "G"},
        {"site": 1, "node": 4, "derived_state": "G"},
        {"site": 2, "node": 5, "derived_state": "G"},
    ],
}

from dnakit.population import count_tree_topologies

result = count_tree_topologies(TABLES, {"A": [0, 1], "B": [2, 3]}, verbose=False)
print(result["tree_index"], result["interval"])
for label, entries in result["counts"].items():
    print(label, [entry["count"] for entry in entries])
```

**示例输出：**

```text
0 [0.0, 100.0]
A [2]
A+B [4]
B [2]
```

## 8) `export_tree_sequence_nexus`

- **作用：** 把树序列写成含 TREES 块的 NEXUS 文件。
- **处理方法：** 调用 `tskit.write_nexus`，`include_trees=True` 时每棵边缘树写成一行 Newick（树名为 `t<左>^<右>`，分类单元为 `n<节点号>`）；`include_alignments=None` 时只要树序列有位点就同时写入样本单倍型序列。已有的 `export_sample_sequences(file_format="nexus")` 固定 `include_trees=False`，只写序列。有多个根的树不能写出。目标文件不得已存在。
- **API：** `export_tree_sequence_nexus(tables[必须], destination[必须], *, include_trees[可选], include_alignments[可选], precision[可选], missing_data_character[可选], max_rows[可选], verbose[可选])`
- **输入：** `tables` 同上；`destination` 是新文件路径；`include_trees`、`include_alignments`（True、False 或 None）、`precision`（1 到 17 或 None）与 `missing_data_character`（单个字符，默认 `?`）。
- **输出：** 字典，含 `destination`、`bytes`、`num_samples`、`num_trees`、`num_sites`、`has_trees_block` 与 `has_data_block`。

**示例代码：**

```python
TABLES = {
    "sequence_length": 100,
    "nodes": [{"time": 0, "sample": True}] * 4 + [{"time": 1}, {"time": 2}, {"time": 3}],
    "edges": [
        {"left": 0, "right": 100, "parent": parent, "child": child}
        for child, parent in [(0, 4), (1, 4), (2, 5), (3, 5), (4, 6), (5, 6)]
    ],
    "sites": [
        {"position": 10, "ancestral_state": "A"},
        {"position": 20, "ancestral_state": "A"},
        {"position": 30, "ancestral_state": "A"},
    ],
    "mutations": [
        {"site": 0, "node": 0, "derived_state": "G"},
        {"site": 1, "node": 4, "derived_state": "G"},
        {"site": 2, "node": 5, "derived_state": "G"},
    ],
}

import tempfile
from pathlib import Path

from dnakit.population import export_tree_sequence_nexus

with tempfile.TemporaryDirectory() as folder:
    target = Path(folder) / "trees.nex"
    result = export_tree_sequence_nexus(TABLES, target, verbose=False)
    lines = target.read_text().splitlines()
print(result["has_trees_block"], result["has_data_block"], result["num_trees"])
print(*lines[-3:], sep="\n")
```

**示例输出：**

```text
True True 1
BEGIN TREES;
  TREE t0^100 = [&R] ((n0:1,n1:1):2,(n2:2,n3:2):1);
END;
```

## 9) `calculate_windowed_divergence`

- **作用：** 按碱基窗口计算两个种群之间每个碱基的核苷酸分歧 dxy。
- **处理方法：** 调用 `allel.windowed_divergence`。已有的 `calculate_sequence_divergence` 对整个区域只给一个值，没有窗口参数。
- **API：** `calculate_windowed_divergence(positions[必须], first_counts[必须], second_counts[必须], *, window_size[必须], start[可选], stop[可选], step[可选], is_accessible[可选], verbose[可选])`
- **输入：** `positions` 是排序的 1 基位置；各个计数数组每个变异一行（每个等位基因一列，至少两列）；`window_size` 是窗口的碱基数；`start`、`stop`、`step` 可选；`is_accessible` 是从位置 1 开始每个碱基一个布尔值的列表，给出时只统计可及碱基。
- **输出：** 字典，含 `window_size`、`step` 与 `windows`；每个窗口含 `start`、`stop`、`sites`、`dxy` 与 `callable_bases`。无法计算的值为 null。

**示例代码：**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_windowed_divergence

result = calculate_windowed_divergence(
    pos, first, second, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], row["callable_bases"], round(row["dxy"], 4))
```

**示例输出：**

```text
1 40 4 40 0.0469
41 80 4 40 0.0344
81 120 4 40 0.0344
```

## 10) `calculate_windowed_watterson_theta`

- **作用：** 按碱基窗口计算每个碱基的 Watterson theta。
- **处理方法：** 调用 `allel.windowed_watterson_theta`。已有的 `calculate_watterson_theta` 对整个序列只给一个值。
- **API：** `calculate_windowed_watterson_theta(positions[必须], allele_counts[必须], *, window_size[必须], start[可选], stop[可选], step[可选], is_accessible[可选], verbose[可选])`
- **输入：** `allele_counts` 是一个种群的计数数组；其余参数同 `calculate_windowed_divergence`。
- **输出：** 字典，含 `window_size`、`step` 与 `windows`；每个窗口含 `start`、`stop`、`sites`、`theta_per_base` 与 `callable_bases`。

**示例代码：**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_windowed_watterson_theta

result = calculate_windowed_watterson_theta(
    pos, first, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], round(row["theta_per_base"], 5))
```

**示例输出：**

```text
1 40 4 0.04091
41 80 4 0.02727
81 120 4 0.02727
```

## 11) `calculate_moving_tajimas_d`

- **作用：** 在连续变异的滑动窗口中计算 Tajima's D。
- **处理方法：** 调用 `allel.moving_tajima_d`。已有的 `calculate_windowed_tajimas_d` 按碱基数划分窗口，窗口的定义不同；窗口内变异少于 `min_sites` 时为 null。
- **API：** `calculate_moving_tajimas_d(allele_counts[必须], *, window_size[必须], start[可选], stop[可选], step[可选], min_sites[可选], verbose[可选])`
- **输入：** 每个变异一行的计数数组；`window_size` 是窗口包含的变异数；`start`、`stop` 是变异下标，`step` 是步长（默认等于 `window_size`）。
- **输出：** 字典，含 `window_size`、`step` 与 `values`（每个窗口一个值）。

**示例代码：**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_moving_tajimas_d

result = calculate_moving_tajimas_d(first, window_size=6, verbose=False)
print(result["window_size"], [round(v, 4) for v in result["values"]])
```

**示例输出：**

```text
6 [0.1677, -0.065]
```

## 12) `calculate_weir_cockerham_fst_per_variant`

- **作用：** 逐变异给出 Weir-Cockerham 方差分量与 Fst。
- **处理方法：** 调用 `allel.weir_cockerham_fst`，把等位基因维度相加：`a` 是种群间分量，`b` 是种群内个体间分量，`c` 是个体内分量，`fst = a / (a + b + c)`。已有的 `calculate_fst(method="weir_cockerham")` 只返回带刀切误差的平均值，没有逐变异输出。
- **API：** `calculate_weir_cockerham_fst_per_variant(genotypes[必须], populations[必须], *, verbose[可选])`
- **输入：** `genotypes` 的形状是（变异、样本、倍性），缺失为 -1；`populations` 是各种群的样本下标列表（至少两个种群，互不重叠）。
- **输出：** 字典，含 `variants`、`populations`（各种群样本数）、逐变异的 `a`、`b`、`c`、`fst`，以及总体的 `fst_ratio_of_sums`。

**示例代码：**

```python
genotypes = [
    [[int((v * (s + 1) + h) % 3 == 0) for h in range(2)] for s in range(4)] for v in range(12)
]
groups = [[0, 1], [2, 3]]

from dnakit.variants import calculate_weir_cockerham_fst_per_variant

result = calculate_weir_cockerham_fst_per_variant(genotypes, groups, verbose=False)
print(result["variants"], result["populations"])
print([round(v, 3) for v in result["fst"][:6]])
print(round(result["fst_ratio_of_sums"], 4))
```

**示例输出：**

```text
12 [2, 2]
[0.0, -0.333, 0.0, 0.0, -0.333, 0.0]
-0.0909
```

## 13) `calculate_windowed_hudson_fst`

- **作用：** 按碱基窗口计算 Hudson Fst（和的比）。
- **处理方法：** 调用 `allel.windowed_hudson_fst`。已有的 `calculate_fst` 没有窗口参数。
- **API：** `calculate_windowed_hudson_fst(positions[必须], first_counts[必须], second_counts[必须], *, window_size[必须], start[可选], stop[可选], step[可选], verbose[可选])`
- **输入：** `positions`、两个种群的计数数组以及窗口参数，含义同 `calculate_windowed_divergence`（没有 `is_accessible`）。
- **输出：** 字典，含 `window_size`、`step` 与 `windows`；每个窗口含 `start`、`stop`、`sites` 与 `fst`。

**示例代码：**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_windowed_hudson_fst

result = calculate_windowed_hudson_fst(
    pos, first, second, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], round(row["fst"], 4))
```

**示例输出：**

```text
1 40 4 0.1111
41 80 4 0.2121
81 120 4 0.0303
```

## 14) `calculate_windowed_weir_cockerham_fst`

- **作用：** 按碱基窗口计算 Weir-Cockerham Fst。
- **处理方法：** 调用 `allel.windowed_weir_cockerham_fst`。已有的 `calculate_fst` 没有窗口参数。
- **API：** `calculate_windowed_weir_cockerham_fst(positions[必须], genotypes[必须], populations[必须], *, window_size[必须], start[可选], stop[可选], step[可选], verbose[可选])`
- **输入：** `positions`、`genotypes` 与 `populations`（同逐变异函数），以及窗口参数。
- **输出：** 字典，含 `window_size`、`step` 与 `windows`；每个窗口含 `start`、`stop`、`sites` 与 `fst`。

**示例代码：**

```python
genotypes = [
    [[int((v * (s + 1) + h) % 3 == 0) for h in range(2)] for s in range(4)] for v in range(12)
]
groups = [[0, 1], [2, 3]]
pos = list(range(10, 130, 10))

from dnakit.variants import calculate_windowed_weir_cockerham_fst

result = calculate_windowed_weir_cockerham_fst(
    pos, genotypes, groups, window_size=40, start=1, stop=120, verbose=False
)
for row in result["windows"]:
    print(row["start"], row["stop"], row["sites"], round(row["fst"], 4))
```

**示例输出：**

```text
1 40 4 -0.0667
41 80 4 -0.1429
81 120 4 -0.0667
```

## 15) `calculate_moving_hudson_fst`

- **作用：** 在连续变异的滑动窗口中计算 Hudson Fst。
- **处理方法：** 调用 `allel.moving_hudson_fst`。已有的 `calculate_fst` 没有窗口参数。
- **API：** `calculate_moving_hudson_fst(first_counts[必须], second_counts[必须], *, window_size[必须], start[可选], stop[可选], step[可选], verbose[可选])`
- **输入：** 每个变异一行的计数数组；`window_size` 是窗口包含的变异数；`start`、`stop` 是变异下标，`step` 是步长（默认等于 `window_size`）。
- **输出：** 字典，含 `window_size`、`step` 与 `values`。

**示例代码：**

```python
pos = list(range(10, 130, 10))
first = [[4 - k, k] for k in [v % 5 for v in range(12)]]
second = [[4 - k, k] for k in [(v * 2 + 1) % 5 for v in range(12)]]

from dnakit.variants import calculate_moving_hudson_fst

result = calculate_moving_hudson_fst(first, second, window_size=4, step=2, verbose=False)
print(result["step"], [round(v, 4) for v in result["values"]])
```

**示例输出：**

```text
2 [0.1111, 0.0667, 0.2121, 0.1282, 0.0303]
```

## 16) `calculate_moving_weir_cockerham_fst`

- **作用：** 在连续变异的滑动窗口中计算 Weir-Cockerham Fst。
- **处理方法：** 调用 `allel.moving_weir_cockerham_fst`。已有的 `calculate_fst` 没有窗口参数。
- **API：** `calculate_moving_weir_cockerham_fst(genotypes[必须], populations[必须], *, window_size[必须], start[可选], stop[可选], step[可选], verbose[可选])`
- **输入：** `genotypes`、`populations` 与滑动窗口参数。
- **输出：** 字典，含 `window_size`、`step` 与 `values`。

**示例代码：**

```python
genotypes = [
    [[int((v * (s + 1) + h) % 3 == 0) for h in range(2)] for s in range(4)] for v in range(12)
]
groups = [[0, 1], [2, 3]]

from dnakit.variants import calculate_moving_weir_cockerham_fst

result = calculate_moving_weir_cockerham_fst(genotypes, groups, window_size=4, step=4, verbose=False)
print(result["step"], [round(v, 4) for v in result["values"]])
```

**示例输出：**

```text
4 [-0.0667, -0.1429, -0.0667]
```

## 17) `calculate_blockwise_patterson_d`

- **作用：** 计算带块刀切标准误与 z 分数的 Patterson's D（ABBA-BABA）。
- **处理方法：** 调用 `allel.blockwise_patterson_d`（`allel.average_patterson_d` 是同一个对象）。四个种群按 A、B、C、D 的顺序给出；每个块含 `block_length` 个连续变异。已有的 `calculate_patterson_d` 只给和的比，没有标准误。
- **API：** `calculate_blockwise_patterson_d(first_counts[必须], second_counts[必须], third_counts[必须], fourth_counts[必须], *, block_length[必须], verbose[可选])`
- **输入：** `first_counts` 到 `fourth_counts` 是四个种群的计数数组（变异数相同）；`block_length` 是块长（变异数）。
- **输出：** 字典，含 `d`、`standard_error`、`z_score`、`block_values`、`jackknife_values`、`variants` 与 `block_length`。

**示例代码：**

```python
def counts(m, c):
    return [[4 - ((v * m + c) % 5), (v * m + c) % 5] for v in range(12)]


a, b, c, d = counts(1, 0), counts(2, 1), counts(1, 2), counts(3, 0)

from dnakit.variants import calculate_blockwise_patterson_d

result = calculate_blockwise_patterson_d(a, b, c, d, block_length=4, verbose=False)
print(round(result["d"], 4), round(result["standard_error"], 4), round(result["z_score"], 4))
print([round(v, 4) for v in result["block_values"]])
```

**示例输出：**

```text
-0.0442 0.1801 -0.2454
[0.0, 0.3404, -0.4364]
```

## 18) `calculate_moving_patterson_d`

- **作用：** 在连续变异的滑动窗口中计算 Patterson's D。
- **处理方法：** 调用 `allel.moving_patterson_d`。已有的 `calculate_patterson_d` 没有窗口参数。
- **API：** `calculate_moving_patterson_d(first_counts[必须], second_counts[必须], third_counts[必须], fourth_counts[必须], *, window_size[必须], start[可选], stop[可选], step[可选], verbose[可选])`
- **输入：** 四个种群的计数数组与滑动窗口参数。
- **输出：** 字典，含 `window_size`、`step` 与 `values`。

**示例代码：**

```python
def counts(m, c):
    return [[4 - ((v * m + c) % 5), (v * m + c) % 5] for v in range(12)]


a, b, c, d = counts(1, 0), counts(2, 1), counts(1, 2), counts(3, 0)

from dnakit.variants import calculate_moving_patterson_d

result = calculate_moving_patterson_d(a, b, c, d, window_size=4, verbose=False)
print(result["window_size"], [round(v, 4) for v in result["values"]])
```

**示例输出：**

```text
4 [0.0, 0.3404, -0.4364]
```

## 19) `calculate_moving_garud_h`

- **作用：** 在连续变异的滑动窗口中计算 Garud 的 H1、H12、H123 与 H2/H1。
- **处理方法：** 调用 `allel.moving_garud_h`，每个窗口独立计算。已有的 `calculate_garud_h` 只对整个数组算一次。
- **API：** `calculate_moving_garud_h(haplotypes[必须], *, window_size[必须], start[可选], stop[可选], step[可选], verbose[可选])`
- **输入：** `haplotypes` 是（变异、单倍型）的 0/1 数组，至少两个单倍型；窗口参数同上。
- **输出：** 字典，含 `window_size`、`step` 以及每个窗口一个值的 `h1`、`h12`、`h123` 与 `h2_h1`。

**示例代码：**

```python
haplotypes = [[int((v * (h + 1)) % 3 == 0) for h in range(6)] for v in range(12)]

from dnakit.variants import calculate_moving_garud_h

result = calculate_moving_garud_h(haplotypes, window_size=6, verbose=False)
print([round(v, 4) for v in result["h1"]])
print([round(v, 4) for v in result["h12"]])
print([round(v, 4) for v in result["h2_h1"]])
```

**示例输出：**

```text
[0.5556, 0.5556]
[1.0, 1.0]
[0.2, 0.2]
```

## 20) `calculate_rogers_huff_r_between`

- **作用：** 计算一组变异与另一组变异之间两两的 Rogers-Huff r。
- **处理方法：** 调用 `allel.rogers_huff_r_between`。已有的 `calculate_rogers_huff_r` 只在同一组变异内部两两配对。每组至多 500 个变异。
- **API：** `calculate_rogers_huff_r_between(alt_counts_a[必须], alt_counts_b[必须], *, verbose[可选])`
- **输入：** `alt_counts_a` 与 `alt_counts_b` 是（变异、样本）的替代等位基因剂量（0、1、2，缺失为 -1），样本必须相同，至少两个样本。
- **输出：** 字典，含 `variants_a`、`variants_b` 与矩阵 `r`（行对应第一组，列对应第二组；常数变异的 r 为 null）。

**示例代码：**

```python
from dnakit.variants import calculate_rogers_huff_r_between

result = calculate_rogers_huff_r_between(
    [[0, 1, 2, 1]], [[0, 1, 2, 1], [2, 1, 0, 1], [1, 1, 1, 1]], verbose=False
)
print(result["variants_a"], result["variants_b"])
print(result["r"])
```

**示例输出：**

```text
1 3
[[1.0, -1.0, None]]
```

## 21) `calculate_randomized_pca`

- **作用：** 对基因型剂量做随机化主成分分析。
- **处理方法：** 去掉单态变异后调用 `allel.randomized_pca`，随机数由 `seed` 固定，同一 `seed` 给出相同结果。随机化求解器近似 `calculate_pca` 的精确 PCA，而 `calculate_pca` 没有求解器参数。`scaler` 可取 `patterson`（默认）、`standard` 或 `none`。
- **API：** `calculate_randomized_pca(alt_counts[必须], *, n_components[可选], iterated_power[可选], seed[可选], scaler[可选], verbose[可选])`
- **输入：** `alt_counts` 是（变异、样本）的剂量 0、1、2，不能有缺失；`n_components` 不超过 50，也不能超过样本数或多态变异数；`iterated_power` 与 `seed` 控制求解器。
- **输出：** 字典，含 `coordinates`（每个样本一行）、`explained_variance_ratio`、`polymorphic_variants`、`n_components`、`seed` 与 `scaler`。

**示例代码：**

```python
dosage = [[(v * (s + 1) + s // 2) % 3 for s in range(6)] for v in range(12)]

from dnakit.variants import calculate_randomized_pca

first = calculate_randomized_pca(dosage, n_components=2, seed=3, verbose=False)
again = calculate_randomized_pca(dosage, n_components=2, seed=3, verbose=False)
print(first["polymorphic_variants"], len(first["coordinates"]), len(first["coordinates"][0]))
print([round(v, 2) for v in first["explained_variance_ratio"]])
print(first == again)
```

**示例输出：**

```text
12 6 2
[0.5, 0.3]
True
```

## 22) `calculate_site_frequency_spectrum_from_counts`

- **作用：** 由双等位基因计数得到位点频谱，可折叠、可缩放。
- **处理方法：** 未折叠时用替代等位基因计数调用 `allel.sfs`（`scale=True` 用 `sfs_scaled`）；`fold=True` 用 `sfs_folded`（`scale=True` 用 `sfs_folded_scaled`）。`spectrum[i]` 是有 `i` 个衍生等位基因的变异数；折叠后有 `sample_size // 2 + 1` 类。已有的 `calculate_site_frequency_spectrum` 读变异记录，没有折叠与缩放选项。
- **API：** `calculate_site_frequency_spectrum_from_counts(allele_counts[必须], *, fold[可选], scale[可选], sample_size[可选], verbose[可选])`
- **输入：** `allele_counts` 的形状是（变异、2），分别为参考与替代（衍生）计数；`sample_size` 默认取各变异总数的最大值，给出时不得小于它。
- **输出：** 字典，含 `variants`、`sample_size`、`folded`、`scaled` 与 `spectrum`（缩放时为浮点数，否则为整数）。

**示例代码：**

```python
from dnakit.variants import calculate_site_frequency_spectrum_from_counts

counts = [[3, 1], [2, 2], [0, 4], [4, 0], [1, 3]]
for fold, scale in [(False, False), (True, False), (False, True), (True, True)]:
    result = calculate_site_frequency_spectrum_from_counts(
        counts, fold=fold, scale=scale, verbose=False
    )
    print(fold, scale, result["sample_size"], result["spectrum"])
```

**示例输出：**

```text
False False 4 [1, 1, 1, 1, 1]
True False 4 [2, 2, 1]
False True 4 [0.0, 1.0, 2.0, 3.0, 4.0]
True True 4 [0.0, 1.5, 1.0]
```

## 23) `calculate_joint_site_frequency_spectrum_from_counts`

- **作用：** 由双等位基因计数得到两个种群的联合位点频谱。
- **处理方法：** 未折叠时调用 `allel.joint_sfs`（`scale=True` 用 `joint_sfs_scaled`）；`fold=True` 调用 `joint_sfs_folded`（`scale=True` 用 `joint_sfs_folded_scaled`），此时每个种群各自折叠到自己的次要等位基因，与 `fold_joint_site_frequency_spectrum` 对合并样本折叠不同。已有的一维频谱包装函数不能给出联合谱。
- **API：** `calculate_joint_site_frequency_spectrum_from_counts(first_counts[必须], second_counts[必须], *, fold[可选], scale[可选], first_sample_size[可选], second_sample_size[可选], verbose[可选])`
- **输入：** `first_counts` 与 `second_counts` 是同一批变异在两个种群中的（变异、2）计数；`first_sample_size` 与 `second_sample_size` 默认取各自总数的最大值。
- **输出：** 字典，含 `variants`、两个样本量、`folded`、`scaled` 与二维的 `spectrum`；`spectrum[i][j]` 对应第一个种群 `i` 个、第二个种群 `j` 个衍生等位基因。

**示例代码：**

```python
from dnakit.variants import calculate_joint_site_frequency_spectrum_from_counts

first = [[3, 1], [2, 2], [0, 4]]
second = [[4, 0], [1, 3], [2, 2]]
plain = calculate_joint_site_frequency_spectrum_from_counts(first, second, verbose=False)
folded = calculate_joint_site_frequency_spectrum_from_counts(first, second, fold=True, verbose=False)
print(plain["first_sample_size"], plain["second_sample_size"])
for row in plain["spectrum"]:
    print(row)
print(folded["spectrum"])
```

**示例输出：**

```text
4 4
[0, 0, 0, 0, 0]
[1, 0, 0, 0, 0]
[0, 0, 0, 1, 0]
[0, 0, 0, 0, 0]
[0, 0, 1, 0, 0]
[[0, 0, 1], [1, 0, 0], [0, 1, 0]]
```

## 24) `phase_trio_step_by_transmission`

- **作用：** 按孟德尔传递对母-父-子三联体做一步定相。
- **处理方法：** `step="progeny"` 调用 `allel.phase_progeny_by_transmission`，在亲本允许时给子代定相；`step="parents"` 调用 `allel.phase_parents_by_transmission`，在 `window_size` 个变异的窗口内用已定相的子代给亲本定相，输入必须是第一步返回的 `genotypes` 与 `is_phased`。已有的 `phase_trio_by_transmission` 总是两步一起运行，不能接受已定相的输入。
- **API：** `phase_trio_step_by_transmission(genotypes[必须], *, step[可选], window_size[可选], is_phased[可选], verbose[可选])`
- **输入：** `genotypes` 的形状是（变异、3、2），样本顺序为母亲、父亲、子代，缺失为 -1；`step` 为 `progeny` 或 `parents`；`parents` 步还需要 `window_size` 与形状为（变异、3）的布尔数组 `is_phased`。
- **输出：** 字典，含 `step`、`variants`、`window_size`、定相后的 `genotypes` 与 `is_phased`。

**示例代码：**

```python
trio = [
    [[(v + s) % 2, (v * (s + 1) + 1) % 2] for s in range(3)] for v in range(12)
]

from dnakit.variants import phase_trio_step_by_transmission

one = phase_trio_step_by_transmission(trio, verbose=False)
two = phase_trio_step_by_transmission(
    one["genotypes"], step="parents", window_size=4, is_phased=one["is_phased"], verbose=False
)
print(sum(row[2] for row in one["is_phased"]), sum(sum(row) for row in one["is_phased"]))
print(sum(sum(row) for row in two["is_phased"]), two["window_size"])
```

**示例输出：**

```text
6 6
18 4
```

## 25) `read_vcf_field_table`

- **作用：** 把本地 VCF 的变异字段读成扁平表格。
- **处理方法：** 调用 `allel.vcf_to_dataframe`（`allel.vcf_to_recarray` 是同一张表的 NumPy 记录数组形式）。`ALT` 与多值 INFO 字段展开成 `alt_number` 列（`ALT_1`、`ALT_2`……）；`fields` 默认是标准固定列，也可写 `variants/DP` 这样的字段名。已有的 `read_vcf` 返回嵌套记录与基因型，没有字段选择与 `alt_number`。
- **API：** `read_vcf_field_table(path[必须], *, fields[可选], alt_number[可选], max_rows[可选], verbose[可选])`
- **输入：** `path` 是本地 VCF；`fields` 是 allel 字段名列表或 None；`alt_number` 在 1 到 20 之间；`max_rows` 是记录数上限（至多 100000）。
- **输出：** 字典，含 `columns`、`num_rows` 与 `rows`（每行一个字典，缺失为 null）。

**示例文件：** `x.vcf`

```text
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	10	rs1	A	G	50	PASS	DP=20;AF=0.5	GT	0/1	1/1
chr1	20	.	C	T,A	.	PASS	DP=7;AF=0.2,0.1	GT	0|1	0/0
chr1	30	rs3	G	C	99	PASS	DP=12	GT	0/0	0/1
```

**示例代码：**

```python
from dnakit.variants import read_vcf_field_table

result = read_vcf_field_table(
    "x.vcf", fields=["CHROM", "POS", "REF", "ALT", "variants/DP"], alt_number=2, verbose=False
)
print(result["columns"], result["num_rows"])
for row in result["rows"]:
    print(row)
```

**示例输出：**

```text
['CHROM', 'POS', 'REF', 'ALT_1', 'ALT_2', 'DP'] 3
{'CHROM': 'chr1', 'POS': 10, 'REF': 'A', 'ALT_1': 'G', 'ALT_2': None, 'DP': 20}
{'CHROM': 'chr1', 'POS': 20, 'REF': 'C', 'ALT_1': 'T', 'ALT_2': 'A', 'DP': 7}
{'CHROM': 'chr1', 'POS': 30, 'REF': 'G', 'ALT_1': 'C', 'ALT_2': None, 'DP': 12}
```

## 26) `vcf_to_csv_table`

- **作用：** 把本地 VCF 的变异字段写成新的 CSV 文件。
- **处理方法：** 调用 `allel.vcf_to_csv`，列与 `read_vcf_field_table` 相同。输出文件必须以 `.csv` 结尾且不得已存在；出错时删除未写完的文件。已有的 `vcf_to_npz_archive` 与 `vcf_to_hdf5_archive` 写数组，不写表格文件。
- **API：** `vcf_to_csv_table(source[必须], output[必须], *, fields[可选], alt_number[可选], verbose[可选])`
- **输入：** `source` 是本地 VCF，`output` 是新 CSV 路径，`fields` 与 `alt_number` 同 `read_vcf_field_table`。
- **输出：** 字典，含 `source`、`output`、`columns`、`num_rows` 与 `bytes`。

**示例代码：**

```python
import tempfile
from pathlib import Path

from dnakit.variants import vcf_to_csv_table

with tempfile.TemporaryDirectory() as folder:
    target = Path(folder) / "table.csv"
    result = vcf_to_csv_table("x.vcf", str(target), fields=["CHROM", "POS", "REF", "ALT"], verbose=False)
    lines = target.read_text().splitlines()
print(result["columns"], result["num_rows"])
print(*lines, sep="\n")
```

**示例输出：**

```text
['CHROM', 'POS', 'REF', 'ALT_1', 'ALT_2', 'ALT_3'] 3
CHROM,POS,REF,ALT_1,ALT_2,ALT_3
chr1,10,A,G,,
chr1,20,C,T,A,
chr1,30,G,C,,
```
