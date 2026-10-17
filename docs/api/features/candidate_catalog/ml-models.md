# 模型训练、评估、解释与导出

本类收录 **48** 项有文档的功能。其中 25 项用普通输入即可计算出结果（`export_torchscript_model` 另需 PyTorch）。5 项需要 DNAKit 之外的东西才能计算，在此之前返回 `status: "not_implemented"`：`attribute`、`mutagenesis`、`explain` 需要在 `options` 中给出 Python 可调用对象（打分函数或模型），`finetune` 与 `generate_sequence_embeddings` 需要本地预训练权重（`checkpoint_path`）或 representation backend。18 项始终返回 `not_implemented`：12 个架构导出（`export_chromatin_model`、`export_cnn_model`、`export_expression_model`、`export_graph_model`、`export_long_context_model`、`export_pangenome_model`、`export_phenotype_model`、`export_protein_model`、`export_rnn_model`、`export_splice_model`、`export_transformer_model`、`export_variant_effect_model`）、两个 ONNX 导出、`pretrain`、`quantize_model`、`run_cross_species_transfer` 和 `serve`。所有函数都通过 DNAKit 命名空间和 Agent/MCP 工具库导出，但 `attribute`、`mutagenesis`、`explain` 与 `export_torchscript_model` 需要 JSON 调用方案无法携带的 Python 可调用对象或 torch 模块，因此只能在 Python 中使用，不能通过 JSON/MCP 使用。下面的示例中，较长的数组和文本用 `...` 缩短。

## 1) attribute

- **作用：** 用用户提供的打分函数，为一条 DNA 序列计算逐位置的计算机诱变（in-silico mutagenesis）归因值。

- **处理方法：** 位置 i 的归因为 `attribution[i] = score(序列) - 位置 i 三种单碱基替换得分的平均值`，得分即所给可调用对象的返回值。DNAKit 不自带模型，打分函数由用户自己提供。未提供打分函数时返回 `not_implemented`。**此函数无法通过 JSON/MCP 使用**：`options['scorer']`（或 `options['model']`）必须是 Python 可调用对象，JSON 调用方案无法携带，只能在 Python 中调用。

- **API：** `dnakit.ml_models.attribute(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：一条无歧义 DNA 序列（A/C/G/T；含歧义碱基会报 `ML_AMBIGUOUS_SEQUENCE`）。`options`：`scorer` 或 `model`（必需，Python 可调用对象，序列 -> 数值）、`positions`（要打分的 0 起始位置，默认全部）、`max_length`（2000）。

- **输出：** 标准封装。`result` 含 `sequence_length`、`method`（`in_silico_mutagenesis`）、`definition`、`reference_score`、`positions` 与 `attributions`（每个位置一个值）。未提供打分函数时：`status: "not_implemented"`，并给出 `required_dependency` 和 `message`。

- **示例代码：**

```python
import json
from dnakit.ml_models import attribute

def a_count(sequence):
    return float(sequence.count("A"))

result = attribute("ACG", options={"scorer": a_count})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "attribute",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 3,
    "sequence_bases": 3,
    "lengths": [3],
    "sha256": "77ef74518ccd6437d9f93ad907bf0bd990f86c2611175f3472aab7852bb1d0e4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_length": 3,
    "method": "in_silico_mutagenesis",
    "definition": "attribution[i] = score(sequence) - mean(score(sequence with base i replaced by each other base))",
    "reference_score": 1.0,
    "positions": [0, 1, 2],
    "attributions": [1.0, -0.3333333333333333, -0.3333333333333333]
  }
}
```

## 2) build

- **作用：** 说明 DNAKit 会为给定模型名与参数创建的（尚未训练的）scikit-learn 估计器。

- **处理方法：** 根据选项构造 `CustomPropertyTrainingConfig`，并通过 DNAKit 的估计器工厂（`dnakit.predictions`）实例化估计器；不拟合、不读取数据。需要 scikit-learn（可用的估计器族如 `random_forest`、`logistic_regression` 等可由 `ml_models.search` 列出）。未知选项名会报 `UNKNOWN_ML_OPTION`。

- **API：** `dnakit.ml_models.build(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`None` 或下列选项组成的映射（与 `options` 合并）。`options`：`task_type`（必需，`classification` 或 `regression`）、`model`（默认 `random_forest`）、`model_params`（估计器参数字典），以及其他 `CustomPropertyTrainingConfig` 字段。

- **输出：** 标准封装。`result` 含 `model`、`task_type`、`estimator_class`（完整类名）、`params`（估计器的完整参数字典）和 `trained: false`。

- **示例代码：**

```python
import json
from dnakit.ml_models import build

result = build(None, options={"task_type": "classification", "model": "random_forest",
                              "model_params": {"n_estimators": 10}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "build",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": "random_forest",
    "task_type": "classification",
    "estimator_class": "sklearn.ensemble._forest.RandomForestClassifier",
    "params": {
      "bootstrap": true,
      "ccp_alpha": 0.0,
      "class_weight": null,
      "criterion": "gini",
      "max_depth": null,
      "max_features": "sqrt",
      "max_leaf_nodes": null,
      "max_samples": null,
      "min_impurity_decrease": 0.0,
      "min_samples_leaf": 1,
      "min_samples_split": 2,
      "min_weight_fraction_leaf": 0.0,
      "monotonic_cst": null,
      "n_estimators": 10,
      "n_jobs": 1,
      "oob_score": false,
      "random_state": 42,
      "verbose": 0,
      "warm_start": false
    },
    "trained": false
  }
}
```

## 3) calculate_prediction_uncertainty

- **作用：** 计算预测不确定性：对类别概率给出熵与间隔；对集成预测给出熵分解或方差。

- **处理方法：** 对每个样本的概率向量计算香农熵 `H = -sum p log p`（默认自然对数）、最小置信度 `1 - max p` 与间隔 `p(第 1) - p(第 2)`。对概率矩阵的集成：预测熵 = 成员平均概率的熵，期望熵 = 各成员熵的平均，互信息 = 预测熵 - 期望熵（认知不确定性部分）。对回归预测的集成：给出成员间的均值、方差、标准差、最小值与最大值（默认总体方差，`ddof=1` 为样本方差）。不做任何校准。

- **API：** `dnakit.ml_models.calculate_prediction_uncertainty(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：含 `probabilities`（样本 x 类别，或成员 x 样本 x 类别）或 `predictions`（成员 x 样本的回归值）之一的映射。`options`：`log_base`（默认 e）、`ddof`（默认 0；取 1 为成员间样本方差）。

- **输出：** 标准封装。`result` 含 `kind`（`classification`、`classification_ensemble` 或 `regression_ensemble`）、`n` 与 `samples`。分类样本含 `predicted_class`、`confidence`、`entropy`、`least_confidence`、`margin`；集成分类样本含 `mean_probabilities`、`predicted_class`、`predictive_entropy`、`expected_entropy`、`mutual_information`；回归集成样本含 `mean`、`variance`、`std`、`min`、`max`。

- **示例代码：**

```python
import json
from dnakit.ml_models import calculate_prediction_uncertainty

result = calculate_prediction_uncertainty({"probabilities": [[0.5, 0.5], [0.9, 0.1]]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "calculate_prediction_uncertainty",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "53d478214e363e76beca7009946ab2ce05820050a42638ac7076d894bf9b271b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "classification",
    "n": 2,
    "samples": [
      {
        "predicted_class": 0,
        "confidence": 0.5,
        "entropy": 0.6931471805599453,
        "least_confidence": 0.5,
        "margin": 0.0
      },
      {
        "predicted_class": 0,
        "confidence": 0.9,
        "entropy": 0.3250829733914482,
        "least_confidence": 0.09999999999999998,
        "margin": 0.8
      }
    ]
  }
}
```

## 4) calibrate

- **作用：** 在二分类得分上拟合概率校准器（Platt 缩放或保序回归），返回校准后的概率与期望校准误差（ECE）。

- **处理方法：** Platt：`p = 1/(1+exp(A*s+B))`，用 Lin-Lin-Weng 平滑目标的牛顿法拟合（测试中与 scikit-learn 的 sigmoid 校准一致）。保序回归：对按升序排列的得分做 pool-adjacent-violators，块之间线性插值，超出拟合范围时截断。ECE 是在 [0, 1] 的 `n_bins` 个等宽分箱上，对 |平均标签 - 平均得分| 按箱内样本数加权的平均；它是在拟合数据上计算的（重代入），因此 `ece_after_on_fit_data` 偏乐观，得分超出 [0, 1] 时 `ece_before` 为 `null`。

- **API：** `dnakit.ml_models.calibrate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：由等长列组成的映射 `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`、此类记录的列表，或 CSV/TSV/JSON/JSONL 路径；列名可用 `true_column`、`pred_column`、`score_column` 修改；`y_score` 必需。`options`：`method`（默认 `platt`，或 `isotonic`）、`positive_label`（默认取较大的标签）、`apply_to`（要变换的得分，默认为拟合所用得分）、`n_bins`（10）。

- **输出：** 标准封装。`result` 含 `method`、`model`（拟合的校准器：`type`、`method`、`classes`、`positive_label`、`parameters`）、`n_fit`、`calibrated`（`apply_to` 对应的概率）、`ece_before`（原始得分的 ECE；只要有任一拟合得分超出 [0, 1] 就为 `null`，因为 ECE 需要概率）、`ece_after_on_fit_data`（始终给出）和 `note`。

- **示例代码：**

```python
import json
from dnakit.ml_models import calibrate

data = {"y_true": [0, 0, 1, 0, 1, 1], "y_score": [0.1, 0.2, 0.35, 0.4, 0.8, 0.9]}
result = calibrate(data, options={"method": "isotonic", "apply_to": [0.15, 0.6]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "calibrate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 74,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e550f9b569b2e26747538240c96c650fe28d96d5680df20f5c22dc7ab0e44a3d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "isotonic",
    "model": {
      "type": "dnakit.calibration",
      "version": 1,
      "method": "isotonic",
      "classes": ["0", "1"],
      "positive_label": "1",
      "parameters": {
        "blocks": [
          {"x_first": 0.1, "x_last": 0.2, "value": 0.0},
          {"x_first": 0.35, "x_last": 0.4, "value": 0.5},
          {"x_first": 0.8, "x_last": 0.9, "value": 1.0}
        ]
      }
    },
    "n_fit": 6,
    "calibrated": [0.0, 0.75],
    "ece_before": 0.27499999999999997,
    "ece_after_on_fit_data": 0.0,
    "note": "ECE is measured on the data used to fit the calibrator (resubstitution); use held-out data for an unbiased estimate."
  }
}
```

## 5) check

- **作用：** 在评估前检查预测列：长度是否一致、有无缺失值与非有限得分、标签是否一致。

- **处理方法：** 读取各列并报告问题而不是抛出异常：长度不一致（`ML_LENGTH_MISMATCH`）与缺失值（`ML_MISSING_VALUES`）会使 `valid` 为假且 `n` 为 null；其余情况下，非有限得分（`NON_FINITE_SCORE`）、回归列中的非数值（`NON_NUMERIC`）和缺失列属于错误，得分超出 [0, 1]、预测标签在 `y_true` 中从未出现属于警告。任务类型按 `auto` 推断：`y_true` 含非整数数值时视为回归。

- **API：** `dnakit.ml_models.check(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：由等长列组成的映射 `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`、此类记录的列表，或 CSV/TSV/JSON/JSONL 路径；列名可用 `true_column`、`pred_column`、`score_column` 修改。`options`：`task`（默认 `auto`，或 `classification`、`regression`）。

- **输出：** 标准封装。`result` 含 `valid`（无错误）、`n`、`task`、`errors`（含 `code`、`message`）与 `warnings`（消息列表）。

- **示例代码：**

```python
import json
from dnakit.ml_models import check

result = check({"y_true": [0, 1], "y_pred": [0, 2], "y_score": [0.2, 1.7]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "check",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 59,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "35d924e53360654aff2b31ce7dc6a74ac892654614a2275375cc9f572ad69679"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": true,
    "n": 2,
    "task": "classification",
    "errors": [],
    "warnings": ["y_score has values outside [0, 1]; they are treated as raw scores, not probabilities.", "y_pred contains labels never seen in y_true: ['2']."]
  }
}
```

## 6) encode

- **作用：** 用 DNAKit 的指纹编码器把 DNA 序列编码为 one-hot 矩阵、整数编码或 k-mer 向量。

- **处理方法：** `one_hot`：每个碱基一行，列为 A、C、G、T。`integer`：A=0、C=1、G=2、T=3（歧义碱基由 `ambiguity_policy` 处理）。`kmer`：全部 4^k 个 k-mer 按 A,C,G,T 字典序的计数、频率或有无，可选规范化（k-mer 与其反向互补合并）。限制：每次调用 `max_sequences` 条序列、每条序列 `max_output_length` 个碱基。

- **API：** `dnakit.ml_models.encode(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：序列字符串、FASTA 文本、序列列表、记录，或表格/FASTA 路径。`options`：`method`（默认 `one_hot`，另有 `integer`、`kmer`）、`k`（3，用于 kmer）、`mode`（`count`、`frequency`、`binary`）、`canonical`（`False`）、`ambiguity_policy`（one_hot：`error` 默认、`zero`、`fractional`；integer：`iupac`、`sentinel`；kmer：`ignore`、`error`）、`sequence_column` / `id_column`（`sequence` / `id`）、`max_sequences`（1000）、`max_output_length`（100000）。

- **输出：** 标准封装。`result` 含 `method`、`n_sequences`、`encodings`（每条序列的 `id`、`length`、`values`），one_hot 与 kmer 另含 `feature_names`。

- **示例代码：**

```python
import json
from dnakit.ml_models import encode

result = encode(["ACGT", "TTA"], options={"method": "integer"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "encode",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 15,
    "sequence_bases": 7,
    "lengths": [4, 3],
    "sha256": "c3ac31c902756f55c7f34b1b8c26595f71d5fe231a9f97f8bce75d5da9d16a84"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "integer",
    "n_sequences": 2,
    "encodings": [
      {"id": "seq1", "length": 4, "values": [0, 1, 2, 3]},
      {"id": "seq2", "length": 3, "values": [3, 3, 0]}
    ]
  }
}
```

## 7) evaluate

- **作用：** 根据预测结果计算分类或回归指标。

- **处理方法：** 分类：准确率、平衡准确率、宏平均/加权的精确率、召回率与 F1、Matthews 相关系数（二分类用公式，多分类用 Gorodkin 的 R_K）、Cohen kappa、混淆矩阵与逐类表；二分类在给出 `y_score` 时另有精确率、召回率/灵敏度、特异度、F1、AUROC（秩统计，并列取平均秩）、AUPRC（平均精确率）与对数损失。回归：MAE、MSE、RMSE、R2、解释方差、Pearson 与 Spearman（无定义的相关系数为 `null`）。任务类型为 `auto`：`y_true` 含非整数数值时视为回归。测试在固定样例上将 AUROC、平均精确率、MCC、F1、kappa、对数损失、R2 等与 scikit-learn 对照。

- **API：** `dnakit.ml_models.evaluate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：由等长列组成的映射 `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`、此类记录的列表，或 CSV/TSV/JSON/JSONL 路径；列名可用 `true_column`、`pred_column`、`score_column` 修改。`options`：`task`（默认 `auto`，或 `classification`、`regression`）、`positive_label`（二分类，默认取较大的标签）、`threshold`（0.5，仅在没有 `y_pred` 而有 `y_score` 时使用）、`classes`（标签顺序）、`metrics`（要返回的指标名子集）。

- **输出：** 标准封装。`result` 含 `task`、`n` 和 `metrics`；分类另含 `per_class`（precision、recall、f1、support）与 `confusion_matrix`（`classes`、`matrix`）。

- **示例代码：**

```python
import json
from dnakit.ml_models import evaluate

data = {"y_true": [0, 0, 1, 1, 1, 0], "y_pred": [0, 1, 1, 1, 0, 0],
        "y_score": [0.1, 0.6, 0.8, 0.7, 0.4, 0.2]}
result = evaluate(data)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "evaluate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 103,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a0420a37d6f513ef723edff4da0739377894cd20eec327bfed038500efc105d7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "task": "classification",
    "n": 6,
    "metrics": {
      "n": 6,
      "accuracy": 0.6666666666666666,
      "balanced_accuracy": 0.6666666666666666,
      "precision_macro": 0.6666666666666666,
      "recall_macro": 0.6666666666666666,
      "f1_macro": 0.6666666666666666,
      "f1_weighted": 0.6666666666666666,
      "mcc": 0.3333333333333334,
      "cohen_kappa": 0.33333333333333326,
      "positive_label": "1",
      "precision": 0.6666666666666666,
      "recall": 0.6666666666666666,
      "sensitivity": 0.6666666666666666,
      "f1": 0.6666666666666666,
      "specificity": 0.6666666666666666,
      "tp": 2,
      "fp": 1,
      "tn": 2,
      "fn": 1,
      "auroc": 0.8888888888888888,
      "auprc": 0.9166666666666666,
      "log_loss": 0.45681733766221466
    },
    "per_class": {
      "0": {
        "precision": 0.6666666666666666,
        "recall": 0.6666666666666666,
        "f1": 0.6666666666666666,
        "support": 3
      },
      "1": {
        "precision": 0.6666666666666666,
        "recall": 0.6666666666666666,
        "f1": 0.6666666666666666,
        "support": 3
      }
    },
    "confusion_matrix": {
      "classes": ["0", "1"],
      "matrix": [
        [2, 1],
        [1, 2]
      ]
    }
  }
}
```

## 8) explain

- **作用：** 用表格数据上的置换特征重要性解释一个模型。

- **处理方法：** 先计算数据上的基线得分（`regression` 用 R2，`classification` 用准确率），然后对每一列做 `n_repeats` 次随机打乱，取得分下降的平均值（`importance_mean` = 基线得分 - 打乱后得分，并给出其标准差）；打乱使用 `random.Random(seed)`。模型不由 DNAKit 提供：`options['model']` 必须是把特征行列表映射为预测的 Python 可调用对象（或带 `predict` 方法的对象），否则返回 `not_implemented`。**此函数无法通过 JSON/MCP 使用**：可调用对象无法放进 JSON 调用方案，只能在 Python 中调用。

- **API：** `dnakit.ml_models.explain(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'X': [[...]], 'y': [...], 'feature_names': [...]}`（`feature_names` 可选，默认 `x0`、`x1` ...）。`options`：`model`（必需，Python 可调用对象或带 `predict` 的对象；`options['scorer']` 可作为 `model` 的别名，两者同时给出时以 `model` 为准）、`task`（默认 `regression`，或 `classification`）、`n_repeats`（5）、`seed`（42）。

- **输出：** 标准封装。`result` 含 `method`（`permutation_importance`）、`metric`（`r2` 或 `accuracy`）、`baseline_score`、`n_repeats` 与 `importances`（每个特征的 `feature`、`importance_mean`、`importance_std`）。未提供模型时：`status: "not_implemented"`，并给出 `required_dependency` 和 `message`。

- **示例代码：**

```python
import json
from dnakit.ml_models import explain

def model(rows):
    return [row[0] for row in rows]      # prediction depends only on feature "f1"

data = {"X": [[1, 5], [2, 3], [3, 8], [4, 1], [5, 9]], "y": [1, 2, 3, 4, 5],
        "feature_names": ["f1", "f2"]}
result = explain(data, options={"model": model, "task": "regression", "n_repeats": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "explain",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9f65a07f1c2350287e76f40e140f8bb4caf3f699469411bb5b635b93a9541e86"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "permutation_importance",
    "metric": "r2",
    "baseline_score": 1.0,
    "n_repeats": 3,
    "importances": [
      {"feature": "f1", "importance_mean": 2.2666666666666666, "importance_std": 0.33993463423951903},
      {"feature": "f2", "importance_mean": 0.0, "importance_std": 0.0}
    ]
  }
}
```

## 9) export_benchmark

- **作用：** 把“模型 x 数据集”的基准测试结果汇总为各数据集内的名次和模型总表，并写出。

- **处理方法：** 在每个数据集内按数值对模型排名（标准竞赛排名：并列者名次相同，下一名次被跳过）。总表为每个模型列出数据集数、平均值、平均名次与获胜次数（名次 1），按平均名次排序。方向由 `higher_is_better` 决定，或按指标名推断（名称含 loss、error、rmse、mae、mse、ece、nll、perplexity 时越小越好）。重复的模型/数据集组合会报 `ML_DUPLICATE_RESULT`。

- **API：** `dnakit.ml_models.export_benchmark(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：长格式记录 `{'model', 'dataset', 'value'}` 或表格路径。`options`：`model_column`、`dataset_column`、`value_column`（默认 `model`、`dataset`、`value`）、`metric`（用于决定方向的指标名）、`higher_is_better`、`format`（默认 `json`：per_dataset + summary；`tsv` / `csv` 写出各数据集的表）、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `higher_is_better`、`format`、`per_dataset`（`dataset`、`model`、`value`、`rank`）、`summary`（`model`、`n_datasets`、`mean_value`、`mean_rank`、`wins`）以及写出文件的字段 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_benchmark

rows = [
    {"model": "A", "dataset": "d1", "value": 0.9}, {"model": "B", "dataset": "d1", "value": 0.8},
    {"model": "A", "dataset": "d2", "value": 0.7}, {"model": "B", "dataset": "d2", "value": 0.8},
]
result = export_benchmark(rows)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_benchmark",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 188,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "57cac65e332b3518094499f796e4b10d0bdc704ae39dab85ee7ff68389cae3ac"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "higher_is_better": true,
    "format": "json",
    "per_dataset": [
      {"dataset": "d1", "model": "A", "value": 0.9, "rank": 1},
      {"dataset": "d1", "model": "B", "value": 0.8, "rank": 2},
      {"dataset": "d2", "model": "A", "value": 0.7, "rank": 2},
      {"dataset": "d2", "model": "B", "value": 0.8, "rank": 1}
    ],
    "summary": [
      {"model": "A", "n_datasets": 2, "mean_value": 0.8, "mean_rank": 1.5, "wins": 1},
      {"model": "B", "n_datasets": 2, "mean_value": 0.8, "mean_rank": 1.5, "wins": 1}
    ],
    "bytes": 655,
    "sha256": "065c29228ad19e6d40ec03896140ed39c702a8c6b6b97418fb8365bec4592984",
    "content": "{\n  \"per_dataset\": [\n    {\n      \"dataset\": \"d1\",\n      \"model\": \"A\",\n      \"value\": 0.9,\n      \"rank\": 1\n    },\n    {\n      \"dataset\": \"d1\",\n      \"model\": \"B\"...",
    "output_path": null
  }
}
```

## 10) export_calibration_model

- **作用：** 在二分类得分上拟合 Platt 或保序校准器，并写成 JSON 模型文件。

- **处理方法：** 拟合方式与 `calibrate` 相同（Platt：对 `p = 1/(1+exp(A*s+B))` 做 Lin-Lin-Weng 平滑目标牛顿拟合；保序：pool-adjacent-violators 分块并线性插值）。JSON 文件含 `type: "dnakit.calibration"`、`version`、`method`、`classes`、`positive_label` 和 `parameters`。

- **API：** `dnakit.ml_models.export_calibration_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：由等长列组成的映射 `{'y_true': [...], 'y_pred': [...], 'y_score': [...]}`、此类记录的列表，或 CSV/TSV/JSON/JSONL 路径；列名可用 `true_column`、`pred_column`、`score_column` 修改；`y_score` 必需。`options`：`method`（`platt` 或 `isotonic`）、`positive_label`、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `model`（校准器）以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_calibration_model

data = {"y_true": [0, 0, 1, 1], "y_score": [0.1, 0.4, 0.6, 0.9]}
result = export_calibration_model(data, options={"method": "isotonic", "output": "calibration.json"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_calibration_model",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 57,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "93a0e7086bac066520556985b76dd7bc916a6e8aba4d95e30579692f1cbcf76b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": {
      "type": "dnakit.calibration",
      "version": 1,
      "method": "isotonic",
      "classes": ["0", "1"],
      "positive_label": "1",
      "parameters": {
        "blocks": [
          {"x_first": 0.1, "x_last": 0.4, "value": 0.0},
          {"x_first": 0.6, "x_last": 0.9, "value": 1.0}
        ]
      }
    },
    "bytes": 351,
    "sha256": "3379ac095b66564832d46e6927ce2c2bd4debb15a3b5cf058557aa74dc2e95bc",
    "output_path": "calibration.json"
  }
}
```

## 11) export_chromatin_model

- **作用：** 把训练好的染色质模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的染色质模型类，也不训练或序列化任何模型；导出染色质模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_chromatin_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_chromatin_model

result = export_chromatin_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_chromatin_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained chromatin model object and its architecture-specific serializer",
  "message": "DNAKit has no chromatin model class to export; nothing was trained or serialised."
}
```

## 12) export_cnn_model

- **作用：** 把训练好的CNN模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的CNN模型类，也不训练或序列化任何模型；导出CNN模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_cnn_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_cnn_model

result = export_cnn_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_cnn_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained cnn model object and its architecture-specific serializer",
  "message": "DNAKit has no cnn model class to export; nothing was trained or serialised."
}
```

## 13) export_dataset

- **作用：** 把记录数据集写成 CSV、TSV、JSON、JSONL 或 FASTA。

- **处理方法：** 读取记录并按所选格式序列化；FASTA 以 `id` 列为序列名、`sequence` 列为序列。文件以原子方式写入（先写临时文件再改名），除非 `overwrite` 为真，否则不会覆盖已有文件。

- **API：** `dnakit.ml_models.export_dataset(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：记录、列映射，或 CSV/TSV/JSON/JSONL/FASTA 路径。`options`：`format`（默认由 `output` 后缀决定，否则 `csv`）、`columns`（列的子集与顺序）、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。FASTA 需要 `sequence` 列。

- **输出：** 标准封装。`result` 含 `format`、`n_records` 以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_dataset

rows = [{"id": "s1", "sequence": "ACGT", "label": 1}, {"id": "s2", "sequence": "GG", "label": 0}]
result = export_dataset(rows, options={"format": "fasta"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_dataset",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 90,
    "sequence_bases": 6,
    "lengths": [4, 2],
    "sha256": "ed916cf32c3cd499c0f54338d0812cffb240170647ba4b99a3b73a7f49bb4bcb"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fasta",
    "n_records": 2,
    "bytes": 16,
    "sha256": "d3dcb573654302517a57f0d7c34c2e0b1706be6d333f0001d397795dd305c803",
    "content": ">s1\nACGT\n>s2\nGG\n",
    "output_path": null
  }
}
```

## 14) export_expression_model

- **作用：** 把训练好的表达模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的表达模型类，也不训练或序列化任何模型；导出表达模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_expression_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_expression_model

result = export_expression_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_expression_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained expression model object and its architecture-specific serializer",
  "message": "DNAKit has no expression model class to export; nothing was trained or serialised."
}
```

## 15) export_graph_model

- **作用：** 把训练好的图模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的图模型类，也不训练或序列化任何模型；导出图模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_graph_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_graph_model

result = export_graph_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_graph_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained graph model object and its architecture-specific serializer",
  "message": "DNAKit has no graph model class to export; nothing was trained or serialised."
}
```

## 16) export_interpretability

- **作用：** 把逐位置归因值（直接提供，或用打分函数通过计算机诱变得到）导出为 TSV 或 JSON。

- **处理方法：** 提供 `attributions` 时按原值写出（`method: "user_supplied"`）；提供一条序列和 `options['scorer']` 时，按 `attribute` 的方式计算归因（得分减去三种替换得分的平均，`method: "in_silico_mutagenesis"`）。只有“直接提供归因值”的形式能通过 JSON/MCP 使用；`scorer` 形式需要 Python 可调用对象，只能在 Python 中使用。

- **API：** `dnakit.ml_models.export_interpretability(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'sequence': 'ACGT', 'attributions': [...]}`，或一条 DNA 序列加 `options['scorer']`（可调用对象，序列 -> 数值）。`options`：`format`（`tsv` 或 `json`）、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `method`、`format`、`n_positions` 以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。TSV 的列为 `position`、`base`、`attribution`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_interpretability

result = export_interpretability({"sequence": "ACG", "attributions": [1.0, -0.25, 0.5]},
                                 options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_interpretability",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 54,
    "sequence_bases": 3,
    "lengths": [3],
    "sha256": "86fc2c5469e64995c2214e7ff45d8407a644e95e93b328c62fdb9630ec18396c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "user_supplied",
    "format": "tsv",
    "n_positions": 3,
    "bytes": 50,
    "sha256": "a76e9bba121ed9b16af1330a282ef9e61756a14d8de97660a794b24e19e03b19",
    "content": "position\tbase\tattribution\n0\tA\t1\n1\tC\t-0.25\n2\tG\t0.5\n",
    "output_path": null
  }
}
```

## 17) export_kmer_embedding

- **作用：** 计算序列的 k-mer 频率（或计数）向量，并写成 TSV 矩阵或 JSON。

- **处理方法：** 每条序列变成长度为 4^k 的向量，k-mer 按 A,C,G,T 顺序（频率 = 计数 / k-mer 总数，或计数，或有无），可选规范化。这是精确的 k-mer 谱，不是学习得到的嵌入。

- **API：** `dnakit.ml_models.export_kmer_embedding(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：与 `encode` 相同（序列字符串、FASTA 文本、列表、记录、表格/FASTA 路径）。`options`：`k`（3）、`mode`（默认 `frequency`，另有 `count`、`binary`）、`canonical`（`False`）、`ambiguity_policy`（`ignore`）、`format`（`tsv` 或 `json`）、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `k`、`mode`、`format`、`dimension`、`n_sequences` 以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。TSV 表头为 `id` 加各 k-mer，每条序列一行。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_kmer_embedding

result = export_kmer_embedding(["AAAA", "ACGT"], options={"k": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_kmer_embedding",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 16,
    "sequence_bases": 8,
    "lengths": [4, 4],
    "sha256": "debb33074aea8f8fcaa086cd144fb49d7cec8ac50769824ba1b3f6924f9c85a6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 1,
    "mode": "frequency",
    "format": "tsv",
    "dimension": 4,
    "n_sequences": 2,
    "bytes": 49,
    "sha256": "1fa9ecba797eaa9bbc97fa4cbb11d67aa7388d29655635d44848b95a5bb31151",
    "content": "id\tA\tC\tG\tT\nseq1\t1\t0\t0\t0\nseq2\t0.25\t0.25\t0.25\t0.25\n",
    "output_path": null
  }
}
```

## 18) export_leaderboard

- **作用：** 按某一指标对模型排名，并把排行榜写成 TSV、CSV、JSON 或 JSONL。

- **处理方法：** 按指标对模型排序，采用标准竞赛排名（数值相同者名次相同，下一名次被跳过）。方向由 `higher_is_better` 指定，或自动判断（指标名含 loss、error、rmse、mae、mse、ece、nll、perplexity 时越小越好）。

- **API：** `dnakit.ml_models.export_leaderboard(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：形如 `[{'model': 'A', 'score': 0.9}, ...]` 的记录，或表格路径。`options`：`metric`（`score`）、`model_column`（`model`）、`higher_is_better`（自动）、`format`（默认 `tsv`，另有 `csv`、`json`、`jsonl`，或由 `output` 后缀决定）、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `metric`、`higher_is_better`、`format`、`n`、`ranking`（`rank`、`model`、指标值）以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_leaderboard

rows = [{"model": "a", "score": 0.9}, {"model": "b", "score": 0.8},
        {"model": "c", "score": 0.8}, {"model": "d", "score": 0.5}]
result = export_leaderboard(rows, options={"format": "json"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_leaderboard",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 4,
    "text_length": 120,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "47f9a912020c2e579d40c9984cb88103254860549918754debc660556057f9db"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "metric": "score",
    "higher_is_better": true,
    "format": "json",
    "n": 4,
    "ranking": [
      {"rank": 1, "model": "a", "score": 0.9},
      {"rank": 2, "model": "b", "score": 0.8},
      {"rank": 2, "model": "c", "score": 0.8},
      {"rank": 4, "model": "d", "score": 0.5}
    ],
    "bytes": 239,
    "sha256": "406c950e5890c02aedb274bc99661459a36734ea4d8a25a6995a96332151e510",
    "content": "[\n  {\n    \"rank\": 1,\n    \"model\": \"a\",\n    \"score\": 0.9\n  },\n  {\n    \"rank\": 2,\n    \"model\": \"b\",\n    \"score\": 0.8\n  },\n  {\n    \"rank\": 2,\n    \"model\": \"c\",\n   ...",
    "output_path": null
  }
}
```

## 19) export_long_context_model

- **作用：** 把训练好的长上下文模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的长上下文模型类，也不训练或序列化任何模型；导出长上下文模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_long_context_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_long_context_model

result = export_long_context_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_long_context_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained long context model object and its architecture-specific serializer",
  "message": "DNAKit has no long context model class to export; nothing was trained or serialised."
}
```

## 20) export_model_onnx

- **作用：** 把训练好的模型导出为 ONNX 格式；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。ONNX 导出需要 `onnx` 与 `onnxruntime` 包（此处未安装）来写出并校验计算图，还需要训练好的模型对象。PyTorch 模块请用 `export_torchscript_model`。此函数与 `export_onnx_model` 是同一功能的重复入口。

- **API：** `dnakit.ml_models.export_model_onnx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_model_onnx

result = export_model_onnx(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_model_onnx",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pip install onnx onnxruntime and a trained model object",
  "message": "ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules."
}
```

## 21) export_onehot_embedding

- **作用：** 对序列做 one-hot 编码，并写成 JSON 或长格式 TSV。

- **处理方法：** 每个碱基变成列为 A、C、G、T 的一行（非 ACGT 碱基按 `ambiguity_policy` 处理：`error` 报错，`zero` 为全零行，`fractional` 把权重平均分给该 IUPAC 代码所代表的碱基，例如 N 各 0.25）。TSV 每个位置一行（`id`、`position`、`A`、`C`、`G`、`T`）；JSON 把每个序列 id 映射到其矩阵。这是固定编码，不是学习得到的嵌入。

- **API：** `dnakit.ml_models.export_onehot_embedding(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：与 `encode` 相同。`options`：`format`（`json` 或 `tsv`，默认由 `output` 后缀决定）、`ambiguity_policy`（默认 `error`，另有 `zero`、`fractional`）、`max_sequences`、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `format`、`n_sequences`、`total_positions` 以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_onehot_embedding

result = export_onehot_embedding([{"id": "s1", "sequence": "AC"}], options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_onehot_embedding",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 32,
    "sequence_bases": 2,
    "lengths": [2],
    "sha256": "6fa4b6fb7d80cb98f276d1517d2e3a9501748dd9111ad544e8a71cc26c91c3f4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "n_sequences": 1,
    "total_positions": 2,
    "bytes": 46,
    "sha256": "5aa43c64f3d9770a1af9bf6a57191c6ad88809a3b0de37e262ebb61c31df0919",
    "content": "id\tposition\tA\tC\tG\tT\ns1\t0\t1\t0\t0\t0\ns1\t1\t0\t1\t0\t0\n",
    "output_path": null
  }
}
```

## 22) export_onnx_model

- **作用：** 把训练好的模型导出为 ONNX 格式；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。ONNX 导出需要 `onnx` 与 `onnxruntime` 包（此处未安装）来写出并校验计算图，还需要训练好的模型对象。PyTorch 模块请用 `export_torchscript_model`。此函数与 `export_model_onnx` 是同一功能的重复入口。

- **API：** `dnakit.ml_models.export_onnx_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_onnx_model

result = export_onnx_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_onnx_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pip install onnx onnxruntime and a trained model object",
  "message": "ONNX export needs the onnx package (not installed) to write and validate the graph; use export_torchscript_model for PyTorch modules."
}
```

## 23) export_pangenome_model

- **作用：** 把训练好的泛基因组模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的泛基因组模型类，也不训练或序列化任何模型；导出泛基因组模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_pangenome_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_pangenome_model

result = export_pangenome_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_pangenome_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained pangenome model object and its architecture-specific serializer",
  "message": "DNAKit has no pangenome model class to export; nothing was trained or serialised."
}
```

## 24) export_phenotype_model

- **作用：** 把训练好的表型模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的表型模型类，也不训练或序列化任何模型；导出表型模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_phenotype_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_phenotype_model

result = export_phenotype_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_phenotype_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained phenotype model object and its architecture-specific serializer",
  "message": "DNAKit has no phenotype model class to export; nothing was trained or serialised."
}
```

## 25) export_protein_model

- **作用：** 把训练好的蛋白质模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的蛋白质模型类，也不训练或序列化任何模型；导出蛋白质模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_protein_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_protein_model

result = export_protein_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_protein_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained protein model object and its architecture-specific serializer",
  "message": "DNAKit has no protein model class to export; nothing was trained or serialised."
}
```

## 26) export_rnn_model

- **作用：** 把训练好的RNN模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的RNN模型类，也不训练或序列化任何模型；导出RNN模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_rnn_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_rnn_model

result = export_rnn_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_rnn_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained rnn model object and its architecture-specific serializer",
  "message": "DNAKit has no rnn model class to export; nothing was trained or serialised."
}
```

## 27) export_splice_model

- **作用：** 把训练好的剪接模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的剪接模型类，也不训练或序列化任何模型；导出剪接模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_splice_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_splice_model

result = export_splice_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_splice_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained splice model object and its architecture-specific serializer",
  "message": "DNAKit has no splice model class to export; nothing was trained or serialised."
}
```

## 28) export_test_split

- **作用：** 写出数据集确定性 train/validation/test 划分中的测试部分。

- **处理方法：** 读取记录并返回 train/validation/test 划分中的一部分。若每条记录都有划分标签（列 `split_column`，默认 `split`；取值 train、validation|valid|val|dev、test，不区分大小写），则直接使用（未知标签报 `ML_UNKNOWN_SPLIT`）。否则用 `random.Random(seed)` 做带种子的随机划分：各部分大小按 `train_ratio`/`validation_ratio`/`test_ratio` 以最大余数法分配（和必须为 1，否则报 `INVALID_SPLIT_RATIOS`）；给出 `group_column` 时同一组的记录留在同一部分。相同输入与种子总得到相同划分，三个 export 合起来恰好覆盖每条记录一次。

- **API：** `dnakit.ml_models.export_test_split(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：记录或表格路径。`options`：`split_column`（`split`）、`train_ratio`（0.8）、`validation_ratio`（0.1）、`test_ratio`（0.1）、`seed`（42）、`group_column`（无）、`format`（默认 `csv`）、`columns`、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `split`、`format`、`n_records`、`counts`（各部分的记录数）、`row_indices`（本部分对应的输入行号，0 起始）以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_test_split

rows = [{"id": f"r{i}", "sequence": "ACGT"} for i in range(10)]
result = export_test_split(rows, options={"seed": 3, "format": "csv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_test_split",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 340,
    "sequence_bases": 40,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "be639db9a82c5bc8d74d72c20d76576f3bc2d901533101ebad790cb7d38ead23"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "split": "test",
    "format": "csv",
    "n_records": 1,
    "counts": {"train": 8, "validation": 1, "test": 1},
    "bytes": 20,
    "sha256": "88bc5ac8ecfcc2658ce1f8c7c983864052ad5cfb27c4ecafbac3d77a5ae4f3cc",
    "content": "id,sequence\nr3,ACGT\n",
    "output_path": null,
    "row_indices": [3]
  }
}
```

## 29) export_tokenizer

- **作用：** 构建 k-mer DNA 分词器词表，并写成 JSON。

- **处理方法：** 词表先列出特殊符号（id 0、1 ...），再按 A,C,G,T 顺序列出全部 4^k 个 k-mer；`k` 最大为 8（更大的值会报配置错误）。JSON 含 `type: "dnakit.kmer_tokenizer"`、模型（`k`、`stride`、`alphabet`、`vocab`）、`special_tokens` 与 `unk_token`。不训练任何分词器，也不使用输入。

- **API：** `dnakit.ml_models.export_tokenizer(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：被忽略（可为 `None`）。`options`：`k`（3，最大 8）、`stride`（1）、`special_tokens`（`['[PAD]', '[UNK]', '[CLS]', '[SEP]', '[MASK]']`）、`output`（JSON 路径；不给则内联返回）、`overwrite`（`False`）。

- **输出：** 标准封装。`result` 含 `k`、`stride`、`vocab_size`、`tokenizer`（词表对象）以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_tokenizer

result = export_tokenizer(None, options={"k": 1})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_tokenizer",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "k": 1,
    "stride": 1,
    "vocab_size": 9,
    "tokenizer": {
      "type": "dnakit.kmer_tokenizer",
      "version": 1,
      "model": {
        "type": "kmer",
        "k": 1,
        "stride": 1,
        "alphabet": "ACGT",
        "vocab": {"[PAD]": 0, "[UNK]": 1, "[CLS]": 2, "[SEP]": 3, "[MASK]": 4, "A": 5, "C": 6, "G": 7, "T": 8}
      },
      "special_tokens": ["[PAD]", "[UNK]", "[CLS]", "[SEP]", "[MASK]"],
      "unk_token": "[UNK]"
    },
    "bytes": 428,
    "sha256": "bbeb600527ff07e8930ef642b5d334cc49fa2ac734568f0c4705e81309079e97",
    "content": "{\n  \"type\": \"dnakit.kmer_tokenizer\",\n  \"version\": 1,\n  \"model\": {\n    \"type\": \"kmer\",\n    \"k\": 1,\n    \"stride\": 1,\n    \"alphabet\": \"ACGT\",\n    \"vocab\": {\n      ...",
    "output_path": null
  }
}
```

## 30) export_torchscript_model

- **作用：** 用示例输入对 PyTorch 模块做追踪（trace），保存为 TorchScript 归档，并验证重新加载后的输出。

- **处理方法：** 在 eval 模式和 `torch.no_grad()` 下执行 `torch.jit.trace`，随后 `save`、`torch.jit.load`，并在示例输入上比较原模型与重新加载后的输出（报告最大绝对差）。模块原先的 training/eval 模式会被恢复。需要 PyTorch（否则返回 `not_implemented`，`required_dependency` 为 `"pip install torch"`）。输入必须是 `torch.nn.Module` 对象，因此此函数无法通过 JSON/MCP 使用，只能在 Python 中调用。ONNX 导出是另一个未实现的函数（`export_onnx_model`）。

- **API：** `dnakit.ml_models.export_torchscript_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`torch.nn.Module`（或 `options['model']`）。`options`：`example_input`（必需，规则的嵌套数字列表，转为 float32 张量）、`output`（必需，`.pt` 路径）、`overwrite`（`False`）。

- **输出：** 标准封装。`result` 含 `output_path`、`bytes`、`input_shape`、`roundtrip_max_abs_difference` 与 `torch_version`。

- **示例代码：**

```python
import json
import torch
from dnakit.ml_models import export_torchscript_model

model = torch.nn.Linear(4, 2)
result = export_torchscript_model(model, options={"example_input": [[0.0, 1.0, 0.0, 0.0]],
                                                  "output": "linear.pt"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_torchscript_model",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 35,
    "lengths": [35],
    "sha256": "6eb7634af5ec8a7ccccd464481ae0dee1f0f1e9237132e1977d422a8a6cc68b9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "output_path": "linear.pt",
    "bytes": 3620,
    "input_shape": [1, 4],
    "roundtrip_max_abs_difference": 0.0,
    "torch_version": "2.6.0+cu118"
  }
}
```

## 31) export_train_split

- **作用：** 写出数据集确定性 train/validation/test 划分中的训练部分。

- **处理方法：** 读取记录并返回 train/validation/test 划分中的一部分。若每条记录都有划分标签（列 `split_column`，默认 `split`；取值 train、validation|valid|val|dev、test，不区分大小写），则直接使用（未知标签报 `ML_UNKNOWN_SPLIT`）。否则用 `random.Random(seed)` 做带种子的随机划分：各部分大小按 `train_ratio`/`validation_ratio`/`test_ratio` 以最大余数法分配（和必须为 1，否则报 `INVALID_SPLIT_RATIOS`）；给出 `group_column` 时同一组的记录留在同一部分。相同输入与种子总得到相同划分，三个 export 合起来恰好覆盖每条记录一次。

- **API：** `dnakit.ml_models.export_train_split(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：记录或表格路径。`options`：`split_column`（`split`）、`train_ratio`（0.8）、`validation_ratio`（0.1）、`test_ratio`（0.1）、`seed`（42）、`group_column`（无）、`format`（默认 `csv`）、`columns`、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `split`、`format`、`n_records`、`counts`（各部分的记录数）、`row_indices`（本部分对应的输入行号，0 起始）以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_train_split

rows = [{"id": f"r{i}", "sequence": "ACGT"} for i in range(10)]
result = export_train_split(rows, options={"seed": 3, "format": "csv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_train_split",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 340,
    "sequence_bases": 40,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "be639db9a82c5bc8d74d72c20d76576f3bc2d901533101ebad790cb7d38ead23"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "split": "train",
    "format": "csv",
    "n_records": 8,
    "counts": {"train": 8, "validation": 1, "test": 1},
    "bytes": 76,
    "sha256": "59391b1a49e78777e0b9a06b303fa4654b760ddbeab878620fd61bbfc93b6766",
    "content": "id,sequence\nr0,ACGT\nr1,ACGT\nr2,ACGT\nr4,ACGT\nr5,ACGT\nr6,ACGT\nr7,ACGT\nr9,ACGT\n",
    "output_path": null,
    "row_indices": [0, 1, 2, 4, 5, 6, "..."]
  }
}
```

## 32) export_transformer_model

- **作用：** 把训练好的Transformer模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的Transformer模型类，也不训练或序列化任何模型；导出Transformer模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_transformer_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_transformer_model

result = export_transformer_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_transformer_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained transformer model object and its architecture-specific serializer",
  "message": "DNAKit has no transformer model class to export; nothing was trained or serialised."
}
```

## 33) export_conformal_uncertainty_model

- **作用：** 在留出的校准集上拟合分割共形（split-conformal）不确定性模型，并写成 JSON。

- **处理方法：** 分割共形预测：非一致性得分为绝对残差 |y - y_hat|（回归）或 1 - p(真实类别)（分类）。保存的分位数是第 k 小的得分，k = ceil((n+1)(1-alpha))；若 k > n（校准集太小），`rank` 与 `quantile` 为 `null`。文件中记录的覆盖率说明是可交换性假设下的边际保证 1 - alpha；DNAKit 不检验可交换性。

- **API：** `dnakit.ml_models.export_conformal_uncertainty_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：`{'y_true': [...], 'y_pred': [...]}`（回归）或 `{'y_true': [...], 'probabilities': [[...]], 'classes': [...]}`（分类）。`options`：`alpha`（0.1）、`output`（JSON 路径；不给则内联返回）、`overwrite`（`False`）。

- **输出：** 标准封装。`result` 含 `model`（`type: "dnakit.split_conformal"`、`task`、`alpha`、`n_calibration`、`rank`、`quantile`、`nonconformity`、`coverage_guarantee`）以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_conformal_uncertainty_model

data = {"y_true": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], "y_pred": [0.0] * 9}
result = export_conformal_uncertainty_model(data, options={"alpha": 0.2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_conformal_uncertainty_model",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 114,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "160c66dde6c7b5f182bfb495e54b0da7b0b0b19f194c5db8f96e4b6428bdb83e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": {
      "type": "dnakit.split_conformal",
      "version": 1,
      "task": "regression",
      "alpha": 0.2,
      "n_calibration": 9,
      "rank": 8,
      "quantile": 8.0,
      "nonconformity": "absolute residual |y - y_hat|",
      "coverage_guarantee": ">= 0.8 marginal coverage under exchangeability"
    },
    "bytes": 275,
    "sha256": "68e78ed9c4245f54b1334992778c392f8bb8643183305896ac3028e9e66eba1b",
    "content": "{\n  \"type\": \"dnakit.split_conformal\",\n  \"version\": 1,\n  \"task\": \"regression\",\n  \"alpha\": 0.2,\n  \"n_calibration\": 9,\n  \"rank\": 8,\n  \"quantile\": 8.0,\n  \"nonconfor...",
    "output_path": null
  }
}
```

## 34) export_validation_split

- **作用：** 写出数据集确定性 train/validation/test 划分中的验证部分。

- **处理方法：** 读取记录并返回 train/validation/test 划分中的一部分。若每条记录都有划分标签（列 `split_column`，默认 `split`；取值 train、validation|valid|val|dev、test，不区分大小写），则直接使用（未知标签报 `ML_UNKNOWN_SPLIT`）。否则用 `random.Random(seed)` 做带种子的随机划分：各部分大小按 `train_ratio`/`validation_ratio`/`test_ratio` 以最大余数法分配（和必须为 1，否则报 `INVALID_SPLIT_RATIOS`）；给出 `group_column` 时同一组的记录留在同一部分。相同输入与种子总得到相同划分，三个 export 合起来恰好覆盖每条记录一次。

- **API：** `dnakit.ml_models.export_validation_split(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：记录或表格路径。`options`：`split_column`（`split`）、`train_ratio`（0.8）、`validation_ratio`（0.1）、`test_ratio`（0.1）、`seed`（42）、`group_column`（无）、`format`（默认 `csv`）、`columns`、`output`（文件路径；不给则不写文件，文本直接返回）、`overwrite`（默认 `False`，已存在的文件除非设为真否则不会被覆盖）。

- **输出：** 标准封装。`result` 含 `split`、`format`、`n_records`、`counts`（各部分的记录数）、`row_indices`（本部分对应的输入行号，0 起始）以及 `bytes`、`sha256`、`output_path`，以及（未给出 `output` 路径时）以 `content` 返回的文本本身。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_validation_split

rows = [{"id": f"r{i}", "sequence": "ACGT"} for i in range(10)]
result = export_validation_split(rows, options={"seed": 3, "format": "csv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_validation_split",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 10,
    "text_length": 340,
    "sequence_bases": 40,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "be639db9a82c5bc8d74d72c20d76576f3bc2d901533101ebad790cb7d38ead23"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "split": "validation",
    "format": "csv",
    "n_records": 1,
    "counts": {"train": 8, "validation": 1, "test": 1},
    "bytes": 20,
    "sha256": "f59306c2734ef4b9cdb82d3f09d0205dad106eae24bd55b63f9d94c0d6815b10",
    "content": "id,sequence\nr8,ACGT\n",
    "output_path": null,
    "row_indices": [8]
  }
}
```

## 35) export_variant_effect_model

- **作用：** 把训练好的变异效应模型导出为可移植文件；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可导出的变异效应模型类，也不训练或序列化任何模型；导出变异效应模型需要其对应架构的序列化器。

- **API：** `dnakit.ml_models.export_variant_effect_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import export_variant_effect_model

result = export_variant_effect_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "export_variant_effect_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained variant effect model object and its architecture-specific serializer",
  "message": "DNAKit has no variant effect model class to export; nothing was trained or serialised."
}
```

## 36) finetune

- **作用：** 在已有本地预训练权重的前提下，通过 `dnakit.predictions.train_custom_property_model` 在带标签文件上微调 DNA 模型（线性探针或全量微调）。

- **处理方法：** 对自定义属性训练器的薄封装，并做参数检查。本调用不会下载权重：若没有 `checkpoint_path`（本地权重）、`representation_backend` 或 `training_encoder`，返回 `not_implemented`。提供后，训练冻结表示上的线性探针（`finetune_method: "linear_probe"`）或微调整个模型（`"full"`），并写出模型目录。其他选项名必须是 `CustomPropertyTrainingConfig` 的字段（否则报 `UNKNOWN_ML_OPTION`）。需要基于 torch 的神经网络运行环境及所选模型自身的依赖。此处未用真实权重运行。

- **API：** `dnakit.ml_models.finetune(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：数据文件路径（CSV/TSV/JSON/FASTA）。`options`：`checkpoint_path`、`representation_backend` 或 `training_encoder`（至少给一个）、`task_type`（必需）、`output_dir`（必需）、`finetune_method`（默认 `linear_probe`，或 `full`）、`model`，以及其他 `CustomPropertyTrainingConfig` 字段。

- **输出：** 没有权重、检查点或编码器时：`status: "not_implemented"`，并给出 `required_dependency`（预训练 DNA 基础模型权重，或 representation backend / training encoder）和 `message`。否则为标准封装，`result.delegate` 与 `result.training`（模型名、目标 schema、划分计数、指标、写出的文件）。

- **示例代码：**

```python
import json
from pathlib import Path
from dnakit.ml_models import finetune

Path("data.csv").write_text("sequence,label\nAAAACC,0\nGGGTTT,1\n")
# no checkpoint_path, representation_backend or training_encoder is given
result = finetune("data.csv", options={"task_type": "classification", "output_dir": "model_out"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "finetune",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 8,
    "sequence_bases": 8,
    "lengths": [8],
    "sha256": "803d8eb0fe6902e3880bee63c537f65f3917f39e743a36703b093591e7054c79"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pretrained DNA foundation-model weights (options['checkpoint_path']) or a representation_backend/training_encoder",
  "message": "Fine-tuning needs a pretrained model; no local checkpoint or encoder was supplied (DNAKit does not download during this call)."
}
```

## 37) generate_sequence_embeddings

- **作用：** 通过 `dnakit.representations.extract_representations`，用 DNAKit 的表示模型为 DNA 序列生成嵌入；需要本地检查点或自行提供的 backend。

- **处理方法：** 委托给表示提取器，对 DNA 语言模型的隐藏状态做池化（`pooling`）。没有 `checkpoint_path`（本地权重）或 `representation_backend`（带 `extract(sequences, show_progress)` 的对象）时返回 `not_implemented`；不下载权重，不联网。它不是 k-mer 谱（k-mer 谱见 `export_kmer_embedding`）。此处未用真实权重运行。

- **API：** `dnakit.ml_models.generate_sequence_embeddings(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：序列字符串、列表、记录或 FASTA/表格路径。`options`：`checkpoint_path` 或 `representation_backend`（至少给一个）、`model`（`ml_models.search` 中的名称）、`model_source_path`、`pooling`、`ambiguity_policy`、`device`、`dtype`、`batch_size`、`max_length`、`max_sequences`（10000）、`sequence_column`、`id_column`。

- **输出：** 没有检查点或 backend 时：`status: "not_implemented"`，并给出 `required_dependency` 和 `message`。否则为标准封装，`result` 含 `delegate`、`model`、`pooling`、`dimension`、`record_ids` 与 `embeddings`（每条序列一个向量）。

- **示例代码：**

```python
import json
from dnakit.ml_models import generate_sequence_embeddings

result = generate_sequence_embeddings(["ACGTACGT", "GGGGCCCC"])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "generate_sequence_embeddings",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 2,
    "text_length": 24,
    "sequence_bases": 16,
    "lengths": [8, 8],
    "sha256": "777448ae73191b73a052271ab7cb66335d35637751e950cb0b24eee1ba9c8c73"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "pretrained DNA foundation-model weights (options['checkpoint_path']) or options['representation_backend']",
  "message": "Sequence embeddings come from a DNA language model; no local checkpoint or backend was supplied."
}
```

## 38) load

- **作用：** 读取由 `ml_models.train` 保存的模型目录的清单（`model.json`），不反序列化任何权重。

- **处理方法：** 解析 `model.json`（至多 10 MB），并检查模型文件（机器学习模型为 `estimator.joblib`，其他为 `head.npz`）是否存在；不会加载 joblib 估计器，因此不会执行任何不受信任的内容。要使用模型，请用 `trust_model=True` 调用 `ml_models.predict`。

- **API：** `dnakit.ml_models.load(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：模型目录路径。没有 `options`。

- **输出：** 标准封装。`result` 含 `model_dir`、`manifest`（任务、schema、特征说明、训练设置、版本）、`artifact`、`artifact_present`、`files`（目录列表）、`weights_loaded: false` 和 `note`。没有 `model.json` 的目录会报 `ML_MODEL_NOT_FOUND`。

- **示例代码：**

```python
import json
from pathlib import Path

Path("data.csv").write_text(
    "sequence,label\nAAAACC,0\nAAACCC,0\nAACCCC,0\nGGGGTT,1\nGGGTTT,1\nGGTTTT,1\nAAAAAA,0\nTTTTTT,1\n")
from dnakit.ml_models import load, train

train("data.csv", options={"task_type": "classification", "output_dir": "model_out",
                           "model": "logistic_regression", "kmer_k": 2, "cv_folds": 2})
result = load("model_out")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "load",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 9,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "74c9778b599cd5ee1b6c3de1b6ea549fe2d68dfe49523bb0e51ad3d99f1100b2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model_dir": "model_out",
    "manifest": {
      "data_mapping": "...",
      "features": "...",
      "finetune_method": "machine_learning",
      "input_dimension": 19,
      "model_name": "logistic_regression",
      "output_dimension": 2,
      "schema_version": "dnakit.custom-property-model.v1",
      "target_schema": "...",
      "task_type": "classification",
      "threshold": 0.5,
      "training": "...",
      "versions": {"joblib": "1.5.3", "numpy": "2.2.6", "python": "3.12.13", "scipy": "1.17.1", "sklearn": "1.8.0"},
      "warnings": []
    },
    "artifact": "estimator.joblib",
    "artifact_present": true,
    "files": ["data_inspection.json", "estimator.joblib", "metrics.json", "model.json", "oof_predictions.jsonl", "splits.csv"],
    "weights_loaded": false,
    "note": "Only the JSON manifest was read; use ml_models.predict with trust_model=True to run the model."
  }
}
```

## 39) mask_dna_tokens

- **作用：** 为掩码语言模型的预处理，随机遮蔽 DNA 序列中一定比例的不重叠 k 碱基片段。

- **处理方法：** 每条序列被切成长度为 `token_size` 的不重叠片段（最后一个片段可能更短）。被遮蔽的片段数为 `round(mask_rate * n_tokens)`，当 `mask_rate > 0` 时至少为 1；遮蔽位置不放回抽样，第 i 条序列的随机数种子为 `seed + i`，因此结果可复现。不涉及任何模型。

- **API：** `dnakit.ml_models.mask_dna_tokens(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：与 `encode` 相同（序列字符串、FASTA 文本、列表、记录、表格/FASTA 路径）。`options`：`token_size`（1）、`mask_rate`（0.15）、`mask_token`（`[MASK]`）、`seed`（42）、`sequence_column`、`id_column`。

- **输出：** 标准封装。`result` 含 `token_size`、`mask_rate`、`mask_token`、`seed` 与 `sequences`；每条序列含 `id`、`tokens`、`masked_tokens`、`masked_positions` 和 `labels`（每个被遮蔽位置的原始片段）。

- **示例代码：**

```python
import json
from dnakit.ml_models import mask_dna_tokens

result = mask_dna_tokens("ACGTACGTAC", options={"token_size": 3, "mask_rate": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "mask_dna_tokens",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 10,
    "lengths": [10],
    "sha256": "6612c0e673bc5b6636f9147556b461347642dcde65d5bc6b39e209d2b6b95af1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "token_size": 3,
    "mask_rate": 0.5,
    "mask_token": "[MASK]",
    "seed": 42,
    "sequences": [
      {
        "id": "seq1",
        "tokens": ["ACG", "TAC", "GTA", "C"],
        "masked_tokens": ["[MASK]", "TAC", "GTA", "[MASK]"],
        "masked_positions": [0, 3],
        "labels": {"0": "ACG", "3": "C"}
      }
    ]
  }
}
```

## 40) mutagenesis

- **作用：** 用用户提供的打分函数，对 DNA 序列的每个单碱基替换打分（计算机饱和诱变）。

- **处理方法：** 对每个选定位置及其他三种碱基，对替换后的序列打分；`delta = 变体得分 - 参考得分`，并按位置给出 delta 的平均值与最大 |delta|。DNAKit 不提供模型，打分函数由用户自己的 Python 函数提供。未提供打分函数时返回 `not_implemented`。**此函数无法通过 JSON/MCP 使用**：`options['scorer']`（或 `options['model']`）必须是 Python 可调用对象，JSON 调用方案无法携带，只能在 Python 中调用。

- **API：** `dnakit.ml_models.mutagenesis(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：一条无歧义 DNA 序列（A/C/G/T；含歧义碱基会报 `ML_AMBIGUOUS_SEQUENCE`）。`options`：`scorer` 或 `model`（必需，Python 可调用对象，序列 -> 数值）、`positions`（0 起始，默认全部）、`max_length`（2000）。

- **输出：** 标准封装。`result` 含 `sequence_length`、`n_variants`、`reference_score`、`effects`（每个变体的 `position`、`ref`、`alt`、`score`、`delta`）与 `per_position`（`position`、`ref`、`mean_delta`、`max_abs_delta`）。未提供打分函数时：`status: "not_implemented"`，并给出 `required_dependency` 和 `message`。

- **示例代码：**

```python
import json
from dnakit.ml_models import mutagenesis

def gc_count(sequence):
    return float(sum(base in "GC" for base in sequence))

result = mutagenesis("AAA", options={"scorer": gc_count, "positions": [1]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "mutagenesis",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 3,
    "sequence_bases": 3,
    "lengths": [3],
    "sha256": "cb1ad2119d8fafb69566510ee712661f9f14b83385006ef92aec47f523a38358"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_length": 3,
    "n_variants": 3,
    "reference_score": 0.0,
    "effects": [
      {"position": 1, "ref": "A", "alt": "C", "score": 1.0, "delta": 1.0},
      {"position": 1, "ref": "A", "alt": "G", "score": 1.0, "delta": 1.0},
      {"position": 1, "ref": "A", "alt": "T", "score": 0.0, "delta": 0.0}
    ],
    "per_position": [
      {"position": 1, "ref": "A", "mean_delta": 0.6666666666666666, "max_abs_delta": 1.0}
    ]
  }
}
```

## 41) predict

- **作用：** 通过 `dnakit.predictions.predict_custom_properties`，用 `ml_models.train` 保存的模型做预测。

- **处理方法：** 加载模型目录，并用训练时相同的特征定义为数据文件中的记录预测（类别标签与类别概率，或回归值）。传统模型以 joblib 文件保存，加载时可能执行代码，因此只有在 `trust_model=True` 时才会加载：只加载你信任的模型。未知选项名会报 `UNKNOWN_ML_OPTION`。

- **API：** `dnakit.ml_models.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：含相同特征列的数据文件路径。`options`：`model_dir`（必需）、`trust_model`（joblib 模型必须为 `True`）、`mapping`（`DataMapping` 字段字典）、`threshold`、`batch_size`、`device`、`output_path`（把预测保存为 csv/tsv/json/parquet）、`overwrite`。

- **输出：** 标准封装。`result` 含 `delegate` 与 `prediction`（`model_dir`、`record_ids`、`label_names`、每条记录预测标签的 `predictions`，分类器另有各类别的 `probabilities`）。

- **示例代码：**

```python
import json
from pathlib import Path

Path("data.csv").write_text(
    "sequence,label\nAAAACC,0\nAAACCC,0\nAACCCC,0\nGGGGTT,1\nGGGTTT,1\nGGTTTT,1\nAAAAAA,0\nTTTTTT,1\n")
from dnakit.ml_models import predict, train

train("data.csv", options={"task_type": "classification", "output_dir": "model_out",
                           "model": "logistic_regression", "kmer_k": 2, "cv_folds": 2})
Path("new.csv").write_text("sequence\nAAAACC\nGGGTTT\n")
result = predict("new.csv", options={"model_dir": "model_out", "trust_model": True})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "predict",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 7,
    "sequence_bases": 6,
    "lengths": [6],
    "sha256": "4d4732a2cb66fff9c0c6d137ffa6452fd04bb9ff620556d7615ee7f746331aa2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "delegate": "dnakit.predictions.predict_custom_properties",
    "prediction": {
      "model_dir": "/work/model_out",
      "record_ids": ["record_1", "record_2"],
      "label_names": ["label"],
      "predictions": [
        {"record_id": "record_1", "label": "0"},
        {"record_id": "record_2", "label": "1"}
      ],
      "probabilities": [
        {
          "record_id": "record_1",
          "label": {"0": 0.9232336713763536, "1": 0.07676632862364648}
        },
        {
          "record_id": "record_2",
          "label": {"0": 0.07620323087640403, "1": 0.923796769123596}
        }
      ]
    }
  }
}
```

## 42) pretrain

- **作用：** 从头预训练 DNA 基础模型；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。自监督预训练需要大规模 DNA 语料、GPU 训练环境和基础模型架构，DNAKit 均不提供。

- **API：** `dnakit.ml_models.pretrain(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import pretrain

result = pretrain(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "pretrain",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "large-scale DNA corpus, GPU training stack and a foundation-model architecture",
  "message": "Self-supervised pre-training of a DNA foundation model is not implemented in DNAKit."
}
```

## 43) quantize_model

- **作用：** 把训练好的神经网络量化为更低精度；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有可量化的神经网络模型对象；需要训练好的神经网络模型和量化后端（`torch.ao` 或 ONNX Runtime）。

- **API：** `dnakit.ml_models.quantize_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import quantize_model

result = quantize_model(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "quantize_model",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained neural model object and a quantisation backend (torch.ao / ONNX Runtime)",
  "message": "Model quantisation is not implemented; DNAKit has no neural model object to quantise."
}
```

## 44) run_cross_species_transfer

- **作用：** 运行跨物种迁移学习实验（在一个物种上训练、在另一个物种上评估）；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。它需要训练好的模型以及带物种标签、迁移协议明确的数据，DNAKit 没有这样的流程。可手动完成：在一个物种的文件上 `train`，再在另一个物种上 `predict` 与 `evaluate`。

- **API：** `dnakit.ml_models.run_cross_species_transfer(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import run_cross_species_transfer

result = run_cross_species_transfer(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "run_cross_species_transfer",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained models and species-labelled data with a defined transfer protocol",
  "message": "Cross-species transfer needs a pretrained model and a species transfer protocol; DNAKit has no such pipeline (use train and predict on the two species files)."
}
```

## 45) search

- **作用：** 搜索 DNAKit 已知的模型目录：传统机器学习估计器、可微调的 DNA 模型、嵌入模型和预训练预测器。

- **处理方法：** 从 `dnakit.predictions` 与 `dnakit.representations` 的注册表构建目录，保留名称、类别或描述中包含查询里每个词（不区分大小写）的条目，可再按 `kind` 与 `task` 过滤。不加载或运行任何模型，也不下载任何内容。

- **API：** `dnakit.ml_models.search(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：查询字符串（空表示全部）。`options`：`query`（可替代 `input`）、`kind`（`machine_learning_estimator`、`dna_model_finetuning`、`dna_embedding_model` 或 `pretrained_prediction_model`）、`task`（`classification`、`regression` 等）、`limit`（50）。

- **输出：** 标准封装。`result` 含 `query`、`kinds`（目录中所有类别）、`n_matches` 与 `matches`（每个条目的 `name`、`kind`、`tasks`、`description`）。

- **示例代码：**

```python
import json
from dnakit.ml_models import search

result = search("random forest")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "search",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 13,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "dc25bdcdacc105874b43917c189e1b76b0dea5a4afe5e9c499e529ea33d48a7d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "query": "random forest",
    "kinds": ["dna_embedding_model", "dna_model_finetuning", "machine_learning_estimator", "pretrained_prediction_model"],
    "n_matches": 1,
    "matches": [
      {
        "name": "random_forest",
        "kind": "machine_learning_estimator",
        "tasks": ["classification", "regression"],
        "description": "scikit-learn estimator"
      }
    ]
  }
}
```

## 46) serve

- **作用：** 通过网络端点提供已训练模型的服务；未实现。

- **处理方法：** DNAKit 未实现此功能：调用不进行任何计算、不写任何文件，始终返回 `status: "not_implemented"`。DNAKit 没有模型服务运行时：调用不会启动服务器、不会打开端口，也不会加载模型。需要 HTTP/gRPC 模型服务器和已部署的模型。

- **API：** `dnakit.ml_models.serve(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 会被接受，仅在 `input_summary` 中汇总；任何 `options` 映射都会被接受并忽略。

- **输出：** 标准封装中的 `status: "not_implemented"`（含 `module`、`operation`、`provenance`、`input_summary`、`warnings`），另有 `required_dependency` 与 `message`，说明缺少什么（见下）。没有 `result`。

- **示例代码：**

```python
import json
from dnakit.ml_models import serve

result = serve(None)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "serve",
  "provenance": {"implementation": "not_implemented", "backend": "local"},
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a model server runtime (HTTP/gRPC) and a deployed model",
  "message": "DNAKit has no model-serving runtime; ml_models.serve does not start a server."
}
```

## 47) train

- **作用：** 通过 `dnakit.predictions.train_custom_property_model`，在带标签的序列/特征文件上训练传统机器学习属性模型（scikit-learn 估计器）。

- **处理方法：** 委托给自定义属性训练器，`finetune_method="machine_learning"`：使用 k-mer（或数值/嵌入）特征，估计器由 `model` 指定（见 `ml_models.search`），可选交叉验证（`cv_folds`）与 train/validation/test 划分，并在全部训练行上重新拟合。写出模型目录（`model.json`、`estimator.joblib`、`metrics.json`、`splits.csv` 等）。这里只训练传统估计器；神经网络微调请用 `ml_models.finetune`。需要 scikit-learn。选项名必须是 `CustomPropertyTrainingConfig` 的字段（否则报 `UNKNOWN_ML_OPTION`）。

- **API：** `dnakit.ml_models.train(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：CSV/TSV/JSON/FASTA 数据文件路径。`options`：`task_type`（必需，`classification` 或 `regression`）、`output_dir`（必需，要写入的模型目录）、`model`（默认 `random_forest`）、`model_params`、`feature_source`（`sequence` k-mer 特征、`numeric`、`embedding`）、`kmer_k`、`mapping`（`DataMapping` 字段字典）、`cv_folds`、`seed`、划分比例、`overwrite`，以及其他 `CustomPropertyTrainingConfig` 字段。

- **输出：** 标准封装。`result` 含 `delegate` 与 `training`（`model_name`、`task_type`、`finetune_method`、`output_dir`、`target_schema`、`split_counts`、含交叉验证各折与折外指标的 `metrics`、`artifacts`、`data_inspection`、`warnings`）。`finetune_method` 不是 `machine_learning` 时返回 `not_implemented`。

- **示例代码：**

```python
import json
from pathlib import Path

Path("data.csv").write_text(
    "sequence,label\nAAAACC,0\nAAACCC,0\nAACCCC,0\nGGGGTT,1\nGGGTTT,1\nGGTTTT,1\nAAAAAA,0\nTTTTTT,1\n")
from dnakit.ml_models import train

result = train("data.csv", options={"task_type": "classification", "output_dir": "model_out",
                           "model": "logistic_regression", "kmer_k": 2, "cv_folds": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "train",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 1,
    "text_length": 8,
    "sequence_bases": 7,
    "lengths": [7],
    "sha256": "803d8eb0fe6902e3880bee63c537f65f3917f39e743a36703b093591e7054c79"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "delegate": "dnakit.predictions.train_custom_property_model",
    "training": {
      "model_name": "logistic_regression",
      "task_type": "classification",
      "finetune_method": "machine_learning",
      "output_dir": "/work/model_out",
      "target_schema": "...",
      "split_counts": {"train": 8, "validation": 0, "test": 0},
      "metrics": {
        "cross_validation": {
          "n_splits": 2,
          "splitter": "StratifiedKFold",
          "seed": 42,
          "folds": "...",
          "out_of_fold": {"accuracy": 1.0, "...": "..."},
          "refit_record_count": 8
        }
      },
      "history": [],
      "artifacts": {
        "manifest": "/work/model_out/model.json",
        "estimator": "/work/model_out/estimator.joblib",
        "metrics": "/work/model_out/metrics.json",
        "data_inspection": "/work/model_out/data_inspection.json",
        "splits": "/work/model_out/splits.csv",
        "oof_predictions": "/work/model_out/oof_predictions.jsonl"
      },
      "data_inspection": "...",
      "warnings": []
    }
  }
}
```

## 48) validate

- **作用：** 检查带标签的 DNA 序列数据集能否用于机器学习：字母表、空字段、重复 id、冲突标签以及划分泄漏。

- **处理方法：** 逐行检查：字符超出 `alphabet`（`INVALID_CHARACTERS`；消息按给定顺序列出字母表，默认为 `ACGTN`）、序列或标签为空、长度超出 `min_length`..`max_length`、id 重复（`DUPLICATE_ID`）、相同序列（不区分大小写）带有不同标签（`CONFLICTING_LABELS`）。给出 `split_column` 时，同一序列出现在不同划分中会报 `SPLIT_LEAKAGE`。标签一致的重复序列只给警告。最多列出 100 条错误。

- **API：** `dnakit.ml_models.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：记录（字典列表）或 CSV/TSV/JSON/JSONL/FASTA 路径。`options`：`sequence_column`（`sequence`）、`label_column`（`label`）、`id_column`（`id`）、`split_column`（无）、`alphabet`（`ACGTN`）、`min_length`（1）、`max_length`（10000000）。

- **输出：** 标准封装。`result` 含 `valid`（无错误）、`errors` 与 `warnings`（含 `row`、`code`、`message`）以及 `summary`（`n_records`、`label_counts`、`length_min`、`length_max`、`duplicate_sequence_rows`）。

- **示例代码：**

```python
import json
from dnakit.ml_models import validate

rows = [
    {"id": "a", "sequence": "ACGT", "label": 1, "split": "train"},
    {"id": "b", "sequence": "ACGT", "label": 1, "split": "test"},
    {"id": "c", "sequence": "AXGT", "label": 0, "split": "train"},
]
result = validate(rows, options={"split_column": "split"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "ml_models",
  "operation": "validate",
  "provenance": {"implementation": "native", "backend": "local"},
  "input_summary": {
    "record_count": 3,
    "text_length": 188,
    "sequence_bases": 12,
    "lengths": [4, 4, 4],
    "sha256": "7a9decb3d90a6d6d42525cd21ab09db13d04d403102ca0dd1daeaf89dfe74bba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "errors": [
      {
        "row": 2,
        "code": "INVALID_CHARACTERS",
        "message": "characters ['X'] are outside the alphabet ACGTN"
      },
      {
        "row": 0,
        "code": "SPLIT_LEAKAGE",
        "message": "identical sequence appears in several splits (rows [0, 1])"
      }
    ],
    "warnings": [
      {"row": null, "code": "DUPLICATE_SEQUENCES", "message": "1 rows repeat an earlier sequence"}
    ],
    "summary": {
      "n_records": 3,
      "label_counts": {"0": 1, "1": 2},
      "length_min": 4,
      "length_max": 4,
      "duplicate_sequence_rows": 1
    }
  }
}
```
