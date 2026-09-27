# 蛋白质特征与结构预测：参数与补充说明

用法和示例见[功能教程](../../api/features/23_prediction_protein_features.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

这些入口调用已有模型或预测服务。模型权重、SDK、服务凭据和使用条件继续由原提供方决定；DNAKit 没有重新训练这些模型。

## 标注蛋白质稀疏特征的生物学类别 {#dnakit.databases.esm_describe_sae_feature-parameters}

`dnakit.databases.esm_describe_sae_feature(*, feature_id, sae_model=None, model=None, n_proteins=None, top_residues_per_protein=None, use_cache=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `feature_id` | integer | 必填 | SAE feature index in &#91;0, 16383&#93;. Use the values returned by ESM_get_sae_features.active_features&#91;&#93;.feature_id. |
| `sae_model` | string | None | SAE checkpoint to label. Cache is keyed on this — different SAEs produce different labels for the same id. {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `model` | string | None | ESMC backbone. {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `n_proteins` | integer | None | Number of panel proteins to run SAE on. Default 10 (the full curated panel). Lower for cheaper exploration at the cost of confidence. {&quot;default&quot;: 10} |
| `top_residues_per_protein` | integer | None | For each protein, take the top-K residues where the target feature activates most strongly, then check their UniProt annotations. Default 3. {&quot;default&quot;: 3} |
| `use_cache` | boolean | None | If true (default) and a cached label exists at ~/.cache/tooluniverse/sae_labels/.../feature_{id}.json, return it instead of rerunning the panel. Set false to force a recompute. {&quot;default&quot;: true} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ESM_API_KEY`。


## 汇总蛋白质指定区间的稀疏特征 {#dnakit.databases.esm_get_region_sae_features-parameters}

`dnakit.databases.esm_get_region_sae_features(*, sequence, start_position, end_position, top_k_features=None, model=None, sae_model=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `start_position` | integer | 必填 | 1-indexed inclusive start of the region of interest. |
| `end_position` | integer | 必填 | 1-indexed inclusive end of the region of interest. Must be &gt;= start_position and within sequence length. Region length cap: 500. |
| `top_k_features` | integer | None | Number of top features (by total &#124;activation&#124; over region) to return. {&quot;default&quot;: 20} |
| `model` | string | None | ESMC base model {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE codebook identifier {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ESM_API_KEY`。


## 提取蛋白质各残基的稀疏特征激活 {#dnakit.databases.esm_get_sae_features-parameters}

`dnakit.databases.esm_get_sae_features(*, sequence, model=None, sae_model=None, position=None, window=None, top_k_per_residue=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Protein amino acid sequence in single-letter code (e.g. &#x27;MEEPQSDPSVEPPLSQETFSDLWKLLPENN&#x27;). Standard 20 amino acids; avoid gaps or non-standard characters. Practical length cap ~2,700 AA. |
| `model` | string | None | ESMC backbone model. Currently SAE is only released for esmc-6b-2024-12 (default). {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE checkpoint name. Default is the layer-60 SAE matching esmc-6b backbone. {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `position` | integer &#124; null | None | Optional 1-indexed residue position. If set, only activations within +/- window residues of this position are returned (use for variant analysis to keep response compact). If null, all residues are returned (large response for long sequences). {&quot;default&quot;: null} |
| `window` | integer | None | Residue window radius around position (only used when position is set). Default 8 = +/- 8 residues = 17-residue window centered on the variant. {&quot;default&quot;: 8} |
| `top_k_per_residue` | integer | None | Cap features returned per residue, sorted by absolute activation. Default 64 = full sparsity (no cap). Set lower (e.g. 16) to reduce response size when no position filter is used. {&quot;default&quot;: 64} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ESM_API_KEY`。


## 计算蛋白质序列的模型对数概率 {#dnakit.databases.esm_score_sequence-parameters}

`dnakit.databases.esm_score_sequence(*, sequence, model=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Protein amino acid sequence in single-letter code to score (e.g. &#x27;MKTAYIAKQRQISFVKSHFSRQ&#x27;). |
| `model` | string | None | ESMC model to use for scoring. {&quot;enum&quot;: &#91;&quot;esmc-300m-2024-12&quot;, &quot;esmc-600m-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-300m-2024-12&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ESM_API_KEY`。


## 使用 AlphaFold2 预测蛋白质三维结构 {#dnakit.databases.nvidia_nim_alphafold2-parameters}

`dnakit.databases.nvidia_nim_alphafold2(*, sequence, algorithm=None, databases=None, e_value=None, iterations=None, relax_prediction=None, skip_template_search=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Amino acid sequence to predict structure for (single letter codes) |
| `algorithm` | string | None | MSA search algorithm. mmseqs2 is faster, jackhmmer is more sensitive {&quot;enum&quot;: &#91;&quot;jackhmmer&quot;, &quot;mmseqs2&quot;&#93;, &quot;default&quot;: &quot;mmseqs2&quot;} |
| `databases` | array | None | Sequence databases for MSA search: uniref90, mgnify, small_bfd {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}, &quot;default&quot;: &#91;&quot;small_bfd&quot;&#93;} |
| `e_value` | number | None | E-value threshold for MSA search {&quot;default&quot;: 0.0001} |
| `iterations` | integer | None | Number of search iterations {&quot;default&quot;: 1} |
| `relax_prediction` | boolean | None | Whether to perform structure relaxation {&quot;default&quot;: false} |
| `skip_template_search` | boolean | None | Skip template-based prediction {&quot;default&quot;: true} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`NVIDIA_API_KEY`。


## 使用 ESMFold 从序列预测蛋白质结构 {#dnakit.databases.nvidia_nim_esmfold-parameters}

`dnakit.databases.nvidia_nim_esmfold(*, sequence, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Amino acid sequence (max 1024 chars, valid: ARNDCQEGHILKMFPSTWYVXBOU) {&quot;maxLength&quot;: 1024} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`NVIDIA_API_KEY`。
