# NUPACK 与 DNA 共享块

## 作用

为 P39/P40 整合批次提供明确的 DNA 操作接口，返回可序列化结果。

## 处理方法

函数调用其对应的上游操作；显式指定本地文件，新输出路径不得已存在。UCSC/黑名单下载仅在主动调用时联网。验证环境缺少 NUPACK 与 BLAST+，相关科学计算流程尚未实测。Enformer 使用与架构和适配器配置一致的本地权重。

## API

`dnakit.secondary_structure.nupack_p39_functions`

| Function | Operation |
| --- | --- |
| `nupack_partition_function` | Complex ensemble partition function and free energy (kcal/mol); Decimal partition values are strings. |
| `nupack_structure_energy` | Free energy in kcal/mol of a specified dot-parens-plus secondary structure. |
| `nupack_sample_structures` | Draw Boltzmann secondary structures using NUPACK's native sampling (num_sample). |
| `nupack_pair_probabilities` | Equilibrium base-pair probability matrix; diagonal entries denote unpaired probabilities. |
| `nupack_mfe_structures` | MFE proxy secondary structures with proxy and stacking-state free energies. |
| `nupack_suboptimal_structures` | Proxy structures within energy_gap kcal/mol of the minimum stacking-state energy. |
| `nupack_ensemble_size` | Number of stacking states, or secondary structures for ensemble='nostacking'. |
| `nupack_design_complex` | Design one complex to a target structure, optionally with IUPAC sequence constraints. |
| `nupack_sequence_distance` | NUPACK sequence distance, including its native run-length notation and nick handling. |
| `nupack_structure_distance` | NUPACK secondary-structure distance, accepting its native structure notations. |
| `nupack_tube_analysis` | Analyze equilibrium complex concentrations of a DNA test tube; optionally save native analysis for reuse. |
| `nupack_complex_concentrations` | Recompute concentrations from a trusted native NUPACK AnalysisResult file without recomputing partition functions. |
| `nupack_design_tubes` | Design multiple target tubes from IUPAC strand constraints, target structures and molar concentrations. |
| `nupack_evaluate_design` | Evaluate fully specified designed DNA strands under the requested physical model and tube targets. |

`dnakit.molecular_lab.geneblocks_p39_functions`

| Function | Operation |
| --- | --- |
| `geneblocks_find_common` | Discover common DNA blocks shared across sequences using upstream all-versus-all BLAST+. |
| `geneblocks_unique_segments` | Coordinates of non-shared segments using CommonBlocks.compute_unique_blocks. |
| `geneblocks_export_csv` | Write the upstream semicolon-separated common-block table to a new file. |
| `geneblocks_common_records` | Export common blocks as JSON representations of upstream sequence records. |
| `geneblocks_unique_records` | Export non-shared segments as upstream sequence records including their DNA. |
| `geneblocks_annotate_common` | Annotate original DNA sequences with discovered common-block feature coordinates. |
| `geneblocks_transfer_annotations` | Transfer strand-aware features between homologous common blocks of annotated DNA records. |
| `geneblocks_plot_common` | Plot shared blocks along DNA sequences; return the image path and underlying block data. |
| `geneblocks_merge_differences` | Merge nearby edit blocks into larger changes using DiffBlocks.merged. |
| `geneblocks_difference_features` | Represent DNA differences as annotated features relative to s1 or s2. |
| `geneblocks_plot_differences` | Plot DNA edit blocks and return the image path and block coordinates. |
| `geneblocks_reconstruct_sequences` | Reconstruct both sequences from a complete user-supplied list of edit and equal blocks. |
| `geneblocks_sequence_record` | Wrap DNA as a sequence record using geneblocks.biotools.sequence_to_record. |
| `geneblocks_circular_equality` | Test whether DNA strings are rotations of one circular sequence using upstream semantics. |
| `geneblocks_extract_subsequence` | Extract a strand-aware 0-based half-open DNA subsequence with geneblocks.Location. |

## 输入

DNA 字符串、JSON 字典/列表、本地基因组文件或模型权重。区间坐标、热力学参数见各函数 docstring；NUPACK 默认材料为 DNA。全部接口支持关键字参数 `verbose`。

## 输出

含完成状态、后端信息和操作结果的字典；文件操作返回输出路径，绘图返回图片路径及对应数据。缺少可选依赖时抛出 `BackendUnavailableError`。各批 manifest 单独记录未验收流程。

**示例代码：**

```python
from dnakit.secondary_structure.nupack_p39_functions import nupack_sequence_distance
from dnakit.exceptions import BackendUnavailableError
try:
    print(nupack_sequence_distance('ACGT', 'ACGA', verbose=False)['distance'])
except BackendUnavailableError:
    print('NUPACK backend unavailable')
```

**示例输出：**

```text
NUPACK backend unavailable
```

示例展示当前环境实际缺少后端时的结果，不代表已经验收 NUPACK 科学计算。接口依据：[NUPACK 官方工具接口](https://docs.nupack.org/utilities/)。

## NUPACK 环自由能（P39 覆盖更正）

**作用：** 计算单个环的自由能和每种堆叠状态的自由能。

**处理方法：** 调用官方明确公开的 `Model.loop_energy` 与 `Model.stack_energies`：https://docs.nupack.org/model/ 。此前将其归为内部助手而跳过的判断有误。

**API：**

```text
nupack_loop_energy(loop, structure, *, model=None, verbose=True)
nupack_loop_stacking_energies(loop, structure, *, model=None, verbose=True)
```

**输入：** 环序列、NUPACK 环结构；`+` 表示链断点。`model` 接受其他包装器相同的物理模型参数，默认 DNA。

**输出：** `free_energy`（kcal/mol），或原生堆叠状态字符串到自由能的 `stacking_energies` 字典。当前验证环境没有 NUPACK，尚未验证数值等价性。

**示例代码：**

```python
from dnakit.secondary_structure.nupack_p39_functions import nupack_loop_energy
from dnakit.exceptions import BackendUnavailableError
try:
    print(nupack_loop_energy('AATT', '....', verbose=False)['free_energy'])
except BackendUnavailableError:
    print('NUPACK unavailable; energy calculation not run.')
```

**示例输出：**

```text
NUPACK unavailable; energy calculation not run.
```
