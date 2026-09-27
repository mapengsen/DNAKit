# 变异效应与蛋白质稳定性：参数与补充说明

用法和示例见[功能教程](../../api/features/23_prediction_variant_effects.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

这些入口调用已有模型或预测服务。模型权重、SDK、服务凭据和使用条件继续由原提供方决定；DNAKit 没有重新训练这些模型。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 评估变异对基因组调控轨道的影响 {#dnakit.databases.alpha_genome_score_variant-parameters}

`dnakit.databases.alpha_genome_score_variant(*, chromosome, position, reference_bases, alternate_bases, output_type=None, organism=None, sequence_length=None, top_n=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chromosome` | string | 必填 | Chromosome, e.g. &#x27;chr22&#x27;. |
| `position` | integer | 必填 | 1-based variant position. |
| `reference_bases` | string | 必填 | Reference allele, e.g. &#x27;A&#x27;. |
| `alternate_bases` | string | 必填 | Alternate allele, e.g. &#x27;C&#x27;. |
| `output_type` | string &#124; null | None | Modality to score: RNA_SEQ (default), ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES, CONTACT_MAPS. |
| `organism` | string &#124; null | None | &#x27;human&#x27; (default) or &#x27;mouse&#x27;. |
| `sequence_length` | string &#124; null | None | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `top_n` | integer &#124; null | None | Number of top &#124;effect&#124; tracks to return (default 20). |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ALPHA_GENOME_API_KEY`。


## 预测氨基酸突变对蛋白质稳定性的影响 {#dnakit.databases.dyna_mut2_predict_stability-parameters}

`dnakit.databases.dyna_mut2_predict_stability(*, operation=None, pdb_id, chain, mutation, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;predict_stability&quot;&#93;, &quot;default&quot;: &quot;predict_stability&quot;} |
| `pdb_id` | string | 必填 | 4-character PDB accession code from RCSB PDB (e.g., &#x27;4HHB&#x27; for hemoglobin, &#x27;1U46&#x27; for ACK1 kinase). The structure will be downloaded automatically. |
| `chain` | string | 必填 | Chain identifier in the PDB structure (e.g., &#x27;A&#x27;). Must be a valid chain in the specified PDB entry. |
| `mutation` | string | 必填 | Single-point mutation in the format WtPosNew where Wt is the wild-type amino acid (1-letter code), Pos is the residue number, and New is the mutant amino acid (1-letter code). Examples: &#x27;V1A&#x27; (Val1Ala), &#x27;E346K&#x27; (Glu346Lys), &#x27;G12D&#x27; (Gly12Asp). |
| `verbose` | bool | True | 运行提示开关。 |


## 使用 ESM-2 评估错义变异效应 {#dnakit.databases.esm2_score_missense_variant-parameters}

`dnakit.databases.esm2_score_missense_variant(*, sequence, position, mutant, wild_type=None, model_id=None, wait_for_model=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Wild-type protein sequence in 1-letter amino-acid code (whitespace is ignored). |
| `position` | integer | 必填 | 1-based residue position of the variant within the sequence. |
| `mutant` | string | 必填 | Mutant amino acid (single 1-letter code, one of the 20 standard residues). |
| `wild_type` | string &#124; null | None | Optional wild-type residue (1-letter). If given, it is validated against the residue at `position` to catch coordinate/sequence mismatches. |
| `model_id` | string &#124; null | None | ESM-2 model to use. Default &#x27;facebook/esm2_t33_650M_UR50D&#x27;. Other sizes: esm2_t30_150M_UR50D, esm2_t12_35M_UR50D, esm2_t6_8M_UR50D. |
| `wait_for_model` | boolean &#124; null | None | If true, block until the model finishes loading on the HF servers instead of returning a retryable &#x27;loading&#x27; status. |
| `verbose` | bool | True | 运行提示开关。 |


## 汇总变异引起的蛋白质特征扰动机制 {#dnakit.databases.esm_explain_variant_mechanism-parameters}

`dnakit.databases.esm_explain_variant_mechanism(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, include_descriptions=None, model=None, sae_model=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `position` | integer | 必填 | 1-indexed residue position of the variant |
| `ref_aa` | string | 必填 | Single-letter wild-type amino acid (must match sequence&#91;position-1&#93;) |
| `alt_aa` | string | 必填 | Single-letter substituted amino acid |
| `window` | integer | None | Residue window centered on the mutation for activation summation. {&quot;default&quot;: 8} |
| `top_k_features` | integer | None | Number of top lost / top gained features to describe and include in the category summary. {&quot;default&quot;: 5} |
| `include_descriptions` | boolean | None | If true (default), call ESM_describe_sae_feature on each top feature to get category labels. If false, return raw feature IDs without categories (cheaper). {&quot;default&quot;: true} |
| `model` | string | None | ESMC base model {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE codebook identifier {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ESM_API_KEY`。


## 批量评估错义变异的稀疏特征扰动 {#dnakit.databases.esm_score_variant_sae_batch-parameters}

`dnakit.databases.esm_score_variant_sae_batch(*, sequence, variants, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `variants` | array | 必填 | List of variants to score. Each variant: {position: 1-indexed int, ref_aa: single letter matching sequence&#91;position-1&#93;, alt_aa: single letter substitution}. {&quot;maxItems&quot;: 100, &quot;items&quot;: {&quot;type&quot;: &quot;object&quot;, &quot;properties&quot;: {&quot;position&quot;: {&quot;type&quot;: &quot;integer&quot;, &quot;description&quot;: &quot;1-indexed residue position&quot;}, &quot;ref_aa&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;description&quot;: &quot;Single-letter wild-type amino acid (must match sequence&#91;position-1&#93;)&quot;}, &quot;alt_aa&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;description&quot;: &quot;Single-letter substituted amino acid&quot;}}, &quot;required&quot;: &#91;&quot;position&quot;, &quot;ref_aa&quot;, &quot;alt_aa&quot;&#93;}} |
| `window` | integer | None | Residue window centered on each variant for activation summation. Default 8 (i.e. positions &#91;pos-8, pos+8&#93;). {&quot;default&quot;: 8} |
| `top_k_features` | integer | None | Number of top features to return per variant (top-K lost + top-K gained). {&quot;default&quot;: 10} |
| `model` | string | None | ESMC base model (must be the 6B model — SAEs are trained against it) {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE codebook identifier {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ESM_API_KEY`。


## 比较参考与突变蛋白质的稀疏特征 {#dnakit.databases.esm_score_variant_sae_disruption-parameters}

`dnakit.databases.esm_score_variant_sae_disruption(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Reference protein sequence (canonical isoform). Single-letter codes, no gaps. The mutant sequence is built internally by substituting alt_aa at position. |
| `position` | integer | 必填 | 1-indexed mutation position. The amino acid at sequence&#91;position-1&#93; must equal ref_aa; otherwise the tool returns an error (the wrong sequence was supplied or the variant notation references a different isoform). |
| `ref_aa` | string | 必填 | Reference amino acid at the mutation position, single-letter code (e.g. &#x27;R&#x27;). Tool validates this matches the supplied sequence. |
| `alt_aa` | string | 必填 | Mutant amino acid, single-letter code (e.g. &#x27;H&#x27; for R175H). |
| `window` | integer | None | Residue window radius around the mutation. Per-feature activations are summed across this window before computing the delta. Default 8 = +/- 8 residues = 17-residue window. {&quot;default&quot;: 8} |
| `top_k_features` | integer | None | Number of top LOST and top GAINED features to return. Default 10. {&quot;default&quot;: 10} |
| `model` | string | None | ESMC backbone (default esmc-6b-2024-12). {&quot;enum&quot;: &#91;&quot;esmc-6b-2024-12&quot;&#93;, &quot;default&quot;: &quot;esmc-6b-2024-12&quot;} |
| `sae_model` | string | None | SAE checkpoint (default layer-60 6B SAE). {&quot;default&quot;: &quot;esmc-6b-2024-12_k64_codebook16384_layer60&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ESM_API_KEY`。


## 使用 Evo 2 评估 DNA 变异效应 {#dnakit.databases.evo2_score_variant-parameters}

`dnakit.databases.evo2_score_variant(*, ref_sequence=None, alt_sequence=None, sequence=None, position=None, reference=None, alternate=None, model=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `ref_sequence` | string &#124; null | None | Reference DNA window (A/C/G/T/N). Use with alt_sequence (same length). |
| `alt_sequence` | string &#124; null | None | Alternate DNA window, same length/centering as ref_sequence. |
| `sequence` | string &#124; null | None | Reference DNA window for point-substitution mode (use with position + alternate). |
| `position` | integer &#124; null | None | 1-based position of the substituted base within `sequence`. |
| `reference` | string &#124; null | None | Optional reference base (single letter) at `position`, validated against `sequence`. |
| `alternate` | string &#124; null | None | Alternate base (single letter) substituted at `position`. |
| `model` | string | None | Hosted Evo 2 model size: evo2-40b (default, most accurate) or evo2-7b (smaller, faster). {&quot;enum&quot;: &#91;&quot;evo2-40b&quot;, &quot;evo2-7b&quot;&#93;, &quot;default&quot;: &quot;evo2-40b&quot;} |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`NVIDIA_API_KEY`。


## 查询变异的最大 SpliceAI 剪接变化分数 {#dnakit.databases.splice_ai_get_max_delta-parameters}

`dnakit.databases.splice_ai_get_max_delta(*, variant, genome=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填 | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) |
| `genome` | string | None | Genome build: 37 or 38 (default: 38) {&quot;enum&quot;: &#91;&quot;37&quot;, &quot;38&quot;&#93;, &quot;default&quot;: &quot;38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 使用 Pangolin 预测变异的剪接影响 {#dnakit.databases.splice_ai_predict_pangolin-parameters}

`dnakit.databases.splice_ai_predict_pangolin(*, variant, genome=None, distance=None, mask=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填 | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) |
| `genome` | string | None | Genome build: 37 or 38 (default: 38) {&quot;enum&quot;: &#91;&quot;37&quot;, &quot;38&quot;&#93;, &quot;default&quot;: &quot;38&quot;} |
| `distance` | integer | None | Distance parameter for model (default: 1000) {&quot;default&quot;: 1000} |
| `mask` | boolean | None | Use masked scores {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |


## 使用 SpliceAI 预测剪接位点变化 {#dnakit.databases.splice_ai_predict_splice-parameters}

`dnakit.databases.splice_ai_predict_splice(*, variant, genome=None, distance=None, mask=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填 | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) or colon-separated |
| `genome` | string | None | Genome build: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 38 {&quot;enum&quot;: &#91;&quot;37&quot;, &quot;38&quot;&#93;, &quot;default&quot;: &quot;38&quot;} |
| `distance` | integer | None | Distance parameter for model (default: 50). Larger values capture more distal effects. {&quot;default&quot;: 50} |
| `mask` | boolean | None | Use masked scores (recommended for variant interpretation). Raw scores better for alternative splicing analysis. {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |
