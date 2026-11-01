# IQ-TREE 3：模型选择、分支支持度、树检验、模拟与定年

`dnakit.phylogeny` 中的这一组函数是对 [IQ-TREE 3](https://github.com/iqtree/iqtree3) 的命令行封装：评估固定树（`--tree-fix`、`--show-lh`）、选择替换模型并合并分区（ModelFinder、`--merge`）、为低分歧基因组推断树（CMAPLE、`--pathogen`）、计算分支支持度（UFBoot、超快 jackknife、bootstrap、jackknife、SH-aLRT、近似贝叶斯、局部 bootstrap 概率、SPRTA）、KH/SH/AU 拓扑检验、似然映射、位点一致性因子、祖先序列重建、对称性检验、比对位点统计、AliSim 模拟、LSD2 定年以及多态感知（PoMo）推断。工作清单中仅方法不同的行合并为带参数的一个函数（`compute_branch_support(methods=...)`、`select_substitution_model(search=..., merge=...)`、`evaluate_tree_likelihood(optimise_branch_lengths=...)`、`compare_tree_topologies(approximately_unbiased=...)`、`compute_site_concordance_factors(method=...)`）。单独的超快 bootstrap 也可通过已有的 `infer_ml_tree` 的 `bootstrap` 参数使用。

测试版本为 IQ-TREE 3.1.4（GPL-2.0，DNAKit 把它作为独立程序启动）。用 `conda install -c bioconda iqtree` 或发布页安装，并把 `iqtree3` 放到 PATH 上，或传入 `iqtree_path`。每个函数都写入一个必须尚不存在的新 `output_directory`，不经过 shell 运行，返回可 JSON 序列化的字典；`outputs` 把前缀之后的文件扩展名（`treefile`、`iqtree`、`contree` 等）映射到路径。值为 `None` 或 `False` 的选项不会传给程序，因此使用 IQ-TREE 的默认值；不给 `model` 时 IQ-TREE 先运行 ModelFinder。示例使用 `verbose=False`。

## 示例数据与示例

**示例文件：** `tree.nwk`

```text
((A:0.1,B:0.12):0.05,((C:0.2,D:0.1):0.07,(E:0.15,F:0.11):0.06):0.04,(G:0.3,H:0.2):0.05);
```

**示例代码：**

```python
from pathlib import Path
import tempfile
from dnakit.phylogeny import iqtree_p11_cli_functions as iq

work = Path(tempfile.mkdtemp())
sim = iq.simulate_alignment(work / "sim", tree="tree.nwk", model="HKY", length=400, seed=7, verbose=False)
print(sim["sequences"], sim["sites"], Path(sim["alignment_files"][0]).name)
stats = iq.alignment_site_statistics(sim["alignment_files"][0], work / "stats", verbose=False)
print(stats["sites"], sorted(stats["counts"]))
```

**示例输出：**

```text
8 400 simulated.fa
400 ['C', 'I', 'U']
```

**示例代码：**

```python
from pathlib import Path
import tempfile
from dnakit.phylogeny import iqtree_p11_cli_functions as iq

work = Path(tempfile.mkdtemp())
sim = iq.simulate_alignment(work / "sim", tree="tree.nwk", model="HKY", length=400, seed=7, verbose=False)
fixed = iq.evaluate_tree_likelihood(sim["alignment_files"][0], "tree.nwk", work / "fit", model="HKY", verbose=False)
print(fixed["version"], round(fixed["log_likelihood"]) < 0)
asr = iq.reconstruct_ancestral_sequences(sim["alignment_files"][0], "tree.nwk", work / "asr", model="HKY", verbose=False)
print(len(asr["nodes"]), asr["sites"])
```

**示例输出：**

```text
3.1.4 True
6 400
```


## 1) `dnakit.phylogeny.iqtree_p11_cli_functions.alignment_site_statistics()`

- **作用：** 把每个比对位点分为恒定、信息位点或无信息位点。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.alignment_site_statistics(alignment: 'str | Path', output_directory: 'str | Path', *, model: 'str' = 'JC', sequence_type: 'str | None' = None, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 2) `dnakit.phylogeny.iqtree_p11_cli_functions.compare_tree_topologies()`

- **作用：** 用 KH、SH 和 AU 检验比较多棵树拓扑。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.compare_tree_topologies(alignment: 'str | Path', trees: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, replicates: 'int' = 10000, approximately_unbiased: 'bool' = True, weighted: 'bool' = False, site_likelihoods: 'bool' = False, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `trees`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 3) `dnakit.phylogeny.iqtree_p11_cli_functions.compute_branch_support()`

- **作用：** 推断最大似然树并用一种或多种方法计算分支支持度（UFBoot、jackknife、bootstrap、SH-aLRT、aBayes、LBP、SPRTA）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.compute_branch_support(alignment: 'str | Path', output_directory: 'str | Path', *, methods: 'Sequence[str]' = ('ufboot',), replicates: 'int | Mapping[str, int]' = 1000, tree: 'str | Path | None' = None, model: 'str | None' = None, jackknife_proportion: 'float | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 4) `dnakit.phylogeny.iqtree_p11_cli_functions.compute_site_concordance_factors()`

- **作用：** 通过采样四分体计算参考树各分支的位点一致性因子（sCF）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.compute_site_concordance_factors(alignment: 'str | Path', tree: 'str | Path', output_directory: 'str | Path', *, quartets: 'int' = 100, method: 'str' = 'likelihood', partitions: 'str | Path | None' = None, model: 'str | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `tree`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 5) `dnakit.phylogeny.iqtree_p11_cli_functions.date_phylogeny()`

- **作用：** 由比对和采样日期用最小二乘法（LSD2）估计时间树和分歧时间。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.date_phylogeny(alignment: 'str | Path', dates: 'str | Path', output_directory: 'str | Path', *, tree: 'str | Path | None' = None, model: 'str | None' = None, tip_date: 'str | None' = None, root_date: 'str | None' = None, confidence_replicates: 'int | None' = None, clock_sd: 'float | None' = None, exclude_outgroup: 'bool' = False, outlier_zscore: 'float | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `dates`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 6) `dnakit.phylogeny.iqtree_p11_cli_functions.evaluate_tree_likelihood()`

- **作用：** 在比对上评估固定树拓扑的对数似然，可选择是否优化枝长。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.evaluate_tree_likelihood(alignment: 'str | Path', tree: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, optimise_branch_lengths: 'bool' = True, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `tree`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 7) `dnakit.phylogeny.iqtree_p11_cli_functions.infer_pathogen_tree()`

- **作用：** 用 CMAPLE 算法为低分歧基因组（病原体/疫情数据）推断树。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.infer_pathogen_tree(alignment: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, force_cmaple: 'bool' = False, alrt_replicates: 'int | None' = None, sprta: 'bool' = False, sprta_zero_length_branches: 'bool' = False, sprta_alternative_placements: 'bool' = False, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 8) `dnakit.phylogeny.iqtree_p11_cli_functions.infer_polymorphism_aware_tree()`

- **作用：** 用多态感知模型（PoMo）由等位基因计数推断群体或物种树。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.infer_polymorphism_aware_tree(counts: 'str | Path', output_directory: 'str | Path', *, substitution_model: 'str' = 'HKY', population_size: 'int | None' = None, sampling: 'str | None' = None, gamma_categories: 'int | None' = None, tree: 'str | Path | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`counts`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 9) `dnakit.phylogeny.iqtree_p11_cli_functions.likelihood_mapping()`

- **作用：** 通过随机四分体的似然映射评估比对的系统发育信号。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.likelihood_mapping(alignment: 'str | Path', output_directory: 'str | Path', *, quartets: 'int' = 10000, clusters: 'str | Path | None' = None, model: 'str | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 10) `dnakit.phylogeny.iqtree_p11_cli_functions.reconstruct_ancestral_sequences()`

- **作用：** 用经验贝叶斯方法重建每个内部节点的祖先序列。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.reconstruct_ancestral_sequences(alignment: 'str | Path', tree: 'str | Path', output_directory: 'str | Path', *, model: 'str | None' = None, min_probability: 'float | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `tree`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 11) `dnakit.phylogeny.iqtree_p11_cli_functions.run_symmetry_tests()`

- **作用：** 用配对对称性检验检查平稳性、可逆性和同质性假设。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.run_symmetry_tests(alignment: 'str | Path', output_directory: 'str | Path', *, only_test: 'bool' = True, remove_partitions: 'str | None' = None, test_type: 'str | None' = None, p_value: 'float | None' = None, keep_zero: 'bool' = False, partitions: 'str | Path | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 12) `dnakit.phylogeny.iqtree_p11_cli_functions.select_substitution_model()`

- **作用：** 选择最优替换模型（ModelFinder），并可合并分区（PartitionFinder 式）。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.select_substitution_model(alignment: 'str | Path', output_directory: 'str | Path', *, search: 'str' = 'MFP', partitions: 'str | Path | None' = None, partition_mode: 'str' = 'proportional', merge: 'bool' = False, merge_algorithm: 'str | None' = None, criterion: 'str | None' = None, model_set: 'str | None' = None, model_subset: 'str | None' = None, rate_models: 'str | None' = None, max_rate_categories: 'int | None' = None, sequence_type: 'str | None' = None, threads: 'int | str' = 1, seed: 'int | None' = None, prefix: 'str' = 'result', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`alignment`, `output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

## 13) `dnakit.phylogeny.iqtree_p11_cli_functions.simulate_alignment()`

- **作用：** 用 AliSim 沿树模拟序列比对。
- **处理方法：** 运行固定的参数列表（不经过 shell）；参数含义、取值范围和工具默认值见英文版同名小节的 Method 说明。
- **API：** `dnakit.phylogeny.iqtree_p11_cli_functions.simulate_alignment(output_directory: 'str | Path', *, tree: 'str | Path | None' = None, random_tree_taxa: 'int | None' = None, random_tree_model: 'str' = 'yh', model: 'str' = 'JC', length: 'int' = 1000, alignments: 'int' = 1, output_format: 'str' = 'fasta', template_alignment: 'str | Path | None' = None, branch_scale: 'float | None' = None, partitions: 'str | Path | None' = None, insertion_deletion_rates: 'Sequence[float] | None' = None, sequence_type: 'str | None' = None, seed: 'int | None' = None, prefix: 'str' = 'simulated', iqtree_path: 'str | None' = None, timeout_seconds: 'float' = 3600.0, verbose: 'bool' = True) -> 'dict[str, Any]'`。
- **输入：** 必填：`output_directory`；其余参数均为仅关键字参数且可选（`None` 或 `False` 表示使用工具自己的默认值）。
- **输出：** 含 `status, backend, operation, version, command, output_directory, outputs, log_tail, provenance` 的字典，另有该函数特有的键（见英文版 Output）。

