# 蛋白质特征与结构预测

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

<!-- BEGIN GENERATED API: model-proteins -->

<span id="model-proteins"></span>

## 1) 标注蛋白质稀疏特征的生物学类别 {#dnakit.databases.esm_describe_sae_feature}

- **作用：** 标注蛋白质稀疏特征的生物学类别。
- **处理方法：** 通过 `ESM` 调用已训练模型，使用所选模型和输入返回预测。
- **API：** `dnakit.databases.esm_describe_sae_feature(*, feature_id, sae_model=None, model=None, n_proteins=None, top_residues_per_protein=None, use_cache=None, verbose=True)`
- **输入：** 必须提供 `feature_id`。 需配置：`ESM_API_KEY`.
- **输出：** 返回字典，包含 `feature_id`, `sae_model`, `category`, `confidence` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.esm_describe_sae_feature(
    feature_id=6404,
    n_proteins=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "feature_id": "<integer>",
  "sae_model": "<string>",
  "category": "<string>",
  "confidence": "<number>",
  "n_proteins_with_activation": "<integer>",
  "n_proteins_analyzed": "<integer>",
  "category_vote_counts": "<object>",
  "supporting_evidence": [
    {
      "protein": "<string>",
      "position_1based": "<integer>",
      "activation": "<number>",
      "uniprot_types": "<array>",
      "uniprot_categories": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_describe_sae_feature-parameters)

## 2) 汇总蛋白质指定区间的稀疏特征 {#dnakit.databases.esm_get_region_sae_features}

- **作用：** 汇总蛋白质指定区间的稀疏特征。
- **处理方法：** 通过 `ESM` 调用已训练模型，使用所选模型和输入返回预测。
- **API：** `dnakit.databases.esm_get_region_sae_features(*, sequence, start_position, end_position, top_k_features=None, model=None, sae_model=None, verbose=True)`
- **输入：** 必须提供 `sequence`, `start_position`, `end_position`。 需配置：`ESM_API_KEY`.
- **输出：** 返回字典，包含 `sequence_length`, `region`, `region_length`, `n_features_active_in_region` 等字段。
- **示例代码：**

```python
import dnakit

example_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'

result = dnakit.databases.esm_get_region_sae_features(
    sequence=example_sequence,
    start_position=10,
    end_position=17,
    top_k_features=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "sequence_length": "<integer>",
  "region": [
    "<value>"
  ],
  "region_length": "<integer>",
  "n_features_active_in_region": "<integer>",
  "top_features": [
    "<value>"
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_get_region_sae_features-parameters)

## 3) 提取蛋白质各残基的稀疏特征激活 {#dnakit.databases.esm_get_sae_features}

- **作用：** 提取蛋白质各残基的稀疏特征激活。
- **处理方法：** 通过 `ESM` 调用已训练模型，使用所选模型和输入返回预测。
- **API：** `dnakit.databases.esm_get_sae_features(*, sequence, model=None, sae_model=None, position=None, window=None, top_k_per_residue=None, verbose=True)`
- **输入：** 必须提供 `sequence`。 需配置：`ESM_API_KEY`.
- **输出：** 返回字典，包含 `sequence_length`, `model`, `sae_model`, `residues_returned` 等字段。
- **示例代码：**

```python
import dnakit

example_sequence = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDR'

result = dnakit.databases.esm_get_sae_features(
    sequence=example_sequence,
    position=175,
    window=8,
    top_k_per_residue=16,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "sequence_length": "<integer>",
  "model": "<string>",
  "sae_model": "<string>",
  "residues_returned": "<integer>",
  "position": "<integer>",
  "window": "<integer>",
  "top_k_per_residue": "<integer>",
  "activations": [
    {
      "residue_idx_1based": "<integer>",
      "active_features": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_get_sae_features-parameters)

## 4) 计算蛋白质序列的模型对数概率 {#dnakit.databases.esm_score_sequence}

- **作用：** 计算蛋白质序列的模型对数概率。
- **处理方法：** 通过 `ESM` 调用已训练模型，使用所选模型和输入返回预测。
- **API：** `dnakit.databases.esm_score_sequence(*, sequence, model=None, verbose=True)`
- **输入：** 必须提供 `sequence`。 需配置：`ESM_API_KEY`.
- **输出：** 返回字典，包含 `status`, `model`, `sequence`, `sequence_length` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.esm_score_sequence(
    sequence='MKTAYIAKQRQISFVKSHFSRQ',
    model='esmc-300m-2024-12',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "model": "<string>",
  "sequence": "<string>",
  "sequence_length": "<integer>",
  "mean_log_probability": "<number>",
  "per_residue_log_probabilities": [
    "<number>"
  ]
}
```

[完整参数与使用说明](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_score_sequence-parameters)

## 5) 使用 AlphaFold2 预测蛋白质三维结构 {#dnakit.databases.nvidia_nim_alphafold2}

- **作用：** 使用 AlphaFold2 预测蛋白质三维结构。
- **处理方法：** 通过 `NvidiaNIM` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.nvidia_nim_alphafold2(*, sequence, algorithm=None, databases=None, e_value=None, iterations=None, relax_prediction=None, skip_template_search=None, verbose=True)`
- **输入：** 必须提供 `sequence`。 需配置：`NVIDIA_API_KEY`.
- **输出：** 返回预测结构的 PDB 坐标及服务提供的置信度和模型指标。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.nvidia_nim_alphafold2(
    sequence='MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回预测结构的 PDB 坐标及服务提供的置信度和模型指标。

[完整参数与使用说明](../../faq/api/23_prediction_protein_features.md#dnakit.databases.nvidia_nim_alphafold2-parameters)

## 6) 使用 ESMFold 从序列预测蛋白质结构 {#dnakit.databases.nvidia_nim_esmfold}

- **作用：** 使用 ESMFold 从序列预测蛋白质结构。
- **处理方法：** 通过 `NvidiaNIM` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.nvidia_nim_esmfold(*, sequence, verbose=True)`
- **输入：** 必须提供 `sequence`。 需配置：`NVIDIA_API_KEY`.
- **输出：** 返回字符串。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.nvidia_nim_esmfold(
    sequence='ARNDCQEGHILKMFPSTWYV',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
"<string>"
```

[完整参数与使用说明](../../faq/api/23_prediction_protein_features.md#dnakit.databases.nvidia_nim_esmfold-parameters)

<!-- END GENERATED API: model-proteins -->
