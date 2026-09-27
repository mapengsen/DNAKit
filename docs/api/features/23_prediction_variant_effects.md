# 变异效应与蛋白质稳定性

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) 长上下文零样本变异效应 {#15}

- **作用：** 计算单 SNV 的 REF、ALT 序列似然及 ALT−REF 变异效应分数。
- **预测方法：** 使用 Evo 2 7B 官方基础模型分别计算完整 REF 和 ALT 上下文的平均序列 log-likelihood，再计算 `alternate − reference`。该值是零样本排序分数，不是临床校准的致病概率。
- **API：** `dnakit.predictions.predict_variant_effects(inputs[必须], config=PropertyPredictionConfig(model="evo2", task="variant_effect")[必须], backend[可选])`
- **输入：** `VariantContext`；REF 与 ALT 必须等长且只有一个 A/C/G/T 位点不同。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="A" * 4096 + "C" + "G" * 4095,
    alternate_sequence="A" * 4096 + "T" + "G" * 4095,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="evo2",
        task="variant_effect",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **示例结果：**

```text
('reference_likelihood', 'alternate_likelihood', 'alternate_minus_reference')
(3,)
```

## 2) 等位基因条件概率变异效应 {#17}

- **作用：** 根据变异位点上游序列，比较下一 token 支持 REF 与 ALT 碱基的概率。
- **预测方法：** 使用 GENERator v2 eukaryote 1.2B 官方权重；最多读取 8,192 bp 上游上下文，去除开头 `N` 并左截断到 6 的倍数，然后计算 `log(p_ref/(p_alt+1e-10))`。该值是零样本排序分数，不是临床校准概率。
- **API：** `dnakit.predictions.predict_variant_effects(inputs[必须], config=PropertyPredictionConfig(model="generator", task="variant_effect", allow_remote_code=True)[必须], backend[可选])`
- **输入：** `VariantContext`；REF 与 ALT 必须等长且只含一个 SNV。模型只使用变异位点之前的参考上下文。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="ACGT" * 1024 + "C" + "G" * 255,
    alternate_sequence="ACGT" * 1024 + "T" + "G" * 255,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="generator",
        task="variant_effect",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **示例结果：**

```text
('reference_probability', 'alternate_probability', 'log_ratio_score')
(3,)
```

## 3) 蛋白稳定性 {#23}

- **作用：** 根据蛋白序列输出官方 ProtStab 任务定义下的稳定性回归值。
- **预测方法：** 调用 LucaOneTasks 官方 `ProtStab` 回归 checkpoint 和 `predict_v1.py` 直接推理，不在 DNAKit 中拟合回归器。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_stability", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="protein"` 的 `BiologicalSequence`。
- **输出：** 返回 `PropertyPredictionResult`；`records` 中每条预测含输出名称 `output_names`、数值 `values` 和相关 metadata。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [
        BiologicalSequence(
            "protein-1",
            "MAMAPRTEINSTRING" * 10,
            "protein",
        )
    ],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_stability",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **示例结果（数值取决于输入和 checkpoint）：**

```text
('prediction',)
(1,)
```

<!-- BEGIN GENERATED API: model-variants -->

<span id="model-variants"></span>

## 4) 评估变异对基因组调控轨道的影响 {#dnakit.databases.alpha_genome_score_variant}

- **作用：** 评估变异对基因组调控轨道的影响。
- **处理方法：** 通过 `AlphaGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.alpha_genome_score_variant(*, chromosome, position, reference_bases, alternate_bases, output_type=None, organism=None, sequence_length=None, top_n=None, verbose=True)`
- **输入：** 必须提供 `chromosome`, `position`, `reference_bases`, `alternate_bases`。 需配置：`ALPHA_GENOME_API_KEY`.
- **输出：** 返回字典，包含 `variant`, `output_type`, `scores` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.alpha_genome_score_variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',
    alternate_bases='C',
    output_type='RNA_SEQ',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": "<string>",
  "output_type": "<string>",
  "scores": [
    "<object>"
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.alpha_genome_score_variant-parameters)

## 5) 预测氨基酸突变对蛋白质稳定性的影响 {#dnakit.databases.dyna_mut2_predict_stability}

- **作用：** 预测氨基酸突变对蛋白质稳定性的影响。
- **处理方法：** 通过 `DynaMut2` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dyna_mut2_predict_stability(*, operation=None, pdb_id, chain, mutation, verbose=True)`
- **输入：** 必须提供 `pdb_id`, `chain`, `mutation`。
- **输出：** 返回字典，包含 `pdb_id`, `chain`, `mutation`, `wild_type` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dyna_mut2_predict_stability(
    operation='predict_stability',
    pdb_id='4HHB',
    chain='A',
    mutation='V1A',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "pdb_id": "<string>",
  "chain": "<string>",
  "mutation": "<string>",
  "wild_type": "<string>",
  "mutant": "<string>",
  "position": "<string>",
  "ddg_prediction": "<number>",
  "effect": "<string>"
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.dyna_mut2_predict_stability-parameters)

## 6) 使用 ESM-2 评估错义变异效应 {#dnakit.databases.esm2_score_missense_variant}

- **作用：** 使用 ESM-2 评估错义变异效应。
- **处理方法：** 通过 `ESM2VariantEffect` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.esm2_score_missense_variant(*, sequence, position, mutant, wild_type=None, model_id=None, wait_for_model=None, verbose=True)`
- **输入：** 必须提供 `sequence`, `position`, `mutant`。
- **输出：** 返回字典，包含 `variant`, `wild_type`, `mutant`, `p_wild_type` 等字段。
- **示例代码：**

```python
import dnakit

example_sequence = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'

result = dnakit.databases.esm2_score_missense_variant(
    sequence=example_sequence,
    position=7,
    wild_type='E',
    mutant='V',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": "<string>",
  "wild_type": "<string>",
  "mutant": "<string>",
  "p_wild_type": "<number>",
  "p_mutant": "<number>",
  "log_likelihood_ratio": "<number>",
  "direction": "<string>"
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm2_score_missense_variant-parameters)

## 7) 汇总变异引起的蛋白质特征扰动机制 {#dnakit.databases.esm_explain_variant_mechanism}

- **作用：** 汇总变异引起的蛋白质特征扰动机制。
- **处理方法：** 通过 `ESM` 调用已训练模型，使用所选模型和输入返回预测。
- **API：** `dnakit.databases.esm_explain_variant_mechanism(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, include_descriptions=None, model=None, sae_model=None, verbose=True)`
- **输入：** 必须提供 `sequence`, `position`, `ref_aa`, `alt_aa`。 需配置：`ESM_API_KEY`.
- **输出：** 返回字典，包含 `variant`, `position`, `window`, `mechanism_summary` 等字段。
- **示例代码：**

```python
import dnakit

example_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'

result = dnakit.databases.esm_explain_variant_mechanism(
    sequence=example_sequence,
    position=12,
    ref_aa='G',
    alt_aa='D',
    top_k_features=3,
    include_descriptions=False,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": "<string>",
  "position": "<integer>",
  "window": "<integer>",
  "mechanism_summary": "<string>",
  "lost_feature_categories": "<object>",
  "gained_feature_categories": "<object>",
  "top_features_lost": [
    "<value>"
  ],
  "top_features_gained": [
    "<value>"
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm_explain_variant_mechanism-parameters)

## 8) 批量评估错义变异的稀疏特征扰动 {#dnakit.databases.esm_score_variant_sae_batch}

- **作用：** 批量评估错义变异的稀疏特征扰动。
- **处理方法：** 通过 `ESM` 调用已训练模型，使用所选模型和输入返回预测。
- **API：** `dnakit.databases.esm_score_variant_sae_batch(*, sequence, variants, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)`
- **输入：** 必须提供 `sequence`, `variants`。 需配置：`ESM_API_KEY`.
- **输出：** 返回字典，包含 `sequence_length`, `n_variants`, `n_succeeded`, `window` 等字段。
- **示例代码：**

```python
import dnakit

example_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'

example_variants = [{'position': 12, 'ref_aa': 'G', 'alt_aa': 'D'},
 {'position': 12, 'ref_aa': 'G', 'alt_aa': 'V'},
 {'position': 11, 'ref_aa': 'A', 'alt_aa': 'T'}]

result = dnakit.databases.esm_score_variant_sae_batch(
    sequence=example_sequence,
    variants=example_variants,
    top_k_features=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "sequence_length": "<integer>",
  "n_variants": "<integer>",
  "n_succeeded": "<integer>",
  "window": "<integer>",
  "results": [
    "<value>"
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm_score_variant_sae_batch-parameters)

## 9) 比较参考与突变蛋白质的稀疏特征 {#dnakit.databases.esm_score_variant_sae_disruption}

- **作用：** 比较参考与突变蛋白质的稀疏特征。
- **处理方法：** 通过 `ESM` 调用已训练模型，使用所选模型和输入返回预测。
- **API：** `dnakit.databases.esm_score_variant_sae_disruption(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)`
- **输入：** 必须提供 `sequence`, `position`, `ref_aa`, `alt_aa`。 需配置：`ESM_API_KEY`.
- **输出：** 返回字典，包含 `variant`, `position`, `window`, `n_unique_features_touched` 等字段。
- **示例代码：**

```python
import dnakit

example_sequence = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDR'

result = dnakit.databases.esm_score_variant_sae_disruption(
    sequence=example_sequence,
    position=175,
    ref_aa='R',
    alt_aa='H',
    window=8,
    top_k_features=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": "<string>",
  "position": "<integer>",
  "window": "<integer>",
  "n_unique_features_touched": "<integer>",
  "top_features_lost": [
    {
      "feature_id": "<integer>",
      "delta": "<number>",
      "ref_activation_sum": "<number>",
      "mut_activation_sum": "<number>"
    }
  ],
  "top_features_gained": [
    {
      "feature_id": "<integer>",
      "delta": "<number>",
      "ref_activation_sum": "<number>",
      "mut_activation_sum": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm_score_variant_sae_disruption-parameters)

## 10) 使用 Evo 2 评估 DNA 变异效应 {#dnakit.databases.evo2_score_variant}

- **作用：** 使用 Evo 2 评估 DNA 变异效应。
- **处理方法：** 通过 `Evo2VariantEffect` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.evo2_score_variant(*, ref_sequence=None, alt_sequence=None, sequence=None, position=None, reference=None, alternate=None, model=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需配置：`NVIDIA_API_KEY`.
- **输出：** 返回字典，包含 `delta_loglik`, `ref_loglik`, `alt_loglik`, `direction` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.evo2_score_variant(
    sequence='ACGTACGTACGTACGTACGTGGGCCCAAATTTACGTACGTACGT',
    position=21,
    reference='G',
    alternate='T',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "delta_loglik": "<number>",
  "ref_loglik": "<number>",
  "alt_loglik": "<number>",
  "direction": "<string>"
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.evo2_score_variant-parameters)

<!-- END GENERATED API: model-variants -->

<!-- BEGIN GENERATED API: variant-splicing -->

<span id="variant-splicing"></span>

## 11) 查询变异的最大 SpliceAI 剪接变化分数 {#dnakit.databases.splice_ai_get_max_delta}

- **作用：** 查询变异的最大 SpliceAI 剪接变化分数。
- **处理方法：** 通过 `SpliceAI` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.splice_ai_get_max_delta(*, variant, genome=None, verbose=True)`
- **输入：** 必须提供 `variant`。
- **输出：** 返回字典，包含 `data`, `source` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.splice_ai_get_max_delta(
    variant='chr8-140300616-T-G',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "variant": "<string>",
    "genome": "<string>",
    "max_delta_score": "<number>",
    "interpretation": "<string>",
    "pathogenicity_threshold": "<string>"
  },
  "source": "<string>"
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.splice_ai_get_max_delta-parameters)

## 12) 使用 Pangolin 预测变异的剪接影响 {#dnakit.databases.splice_ai_predict_pangolin}

- **作用：** 使用 Pangolin 预测变异的剪接影响。
- **处理方法：** 通过 `SpliceAI` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.splice_ai_predict_pangolin(*, variant, genome=None, distance=None, mask=None, verbose=True)`
- **输入：** 必须提供 `variant`。
- **输出：** 返回字典，包含 `data`, `source` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.splice_ai_predict_pangolin(
    variant='chr8-140300616-T-G',
    genome='38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "variant": "<string>",
    "genome": "<string>",
    "scores": [
      "<value>"
    ]
  },
  "source": "<string>"
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.splice_ai_predict_pangolin-parameters)

## 13) 使用 SpliceAI 预测剪接位点变化 {#dnakit.databases.splice_ai_predict_splice}

- **作用：** 使用 SpliceAI 预测剪接位点变化。
- **处理方法：** 通过 `SpliceAI` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.splice_ai_predict_splice(*, variant, genome=None, distance=None, mask=None, verbose=True)`
- **输入：** 必须提供 `variant`。
- **输出：** 返回字典，包含 `data`, `source` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.splice_ai_predict_splice(
    variant='chr8-140300616-T-G',
    genome='38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "variant": "<string>",
    "genome": "<string>",
    "scores": [
      "<value>"
    ],
    "max_delta_score": "<number>",
    "interpretation": "<string>"
  },
  "source": "<string>"
}
```

[完整参数与使用说明](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.splice_ai_predict_splice-parameters)

<!-- END GENERATED API: variant-splicing -->
