# 自定义模型训练与预测

依赖、权重和运行配置见[FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 1) 机器学习分类与回归 {#_2}

- **作用：** 用带标签的数据文件训练分类或回归模型，并预测新样本。
- **训练方法：** 提取序列特征并拟合随机森林。
- **API：** `dnakit.predictions.train_custom_property_model(source, config=...)`; `dnakit.predictions.predict_custom_properties(source, model_dir=...)`
- **输入：** CSV 的 `sequence`、`target` 列；示例为回归，需先安装 `.[ml]`。仅为自己训练或可信来源的模型设置 `trust_model=True`。
- **输出：** 返回 `CustomPropertyTrainingResult`，包含训练产物路径和评估指标。
- **示例代码：**

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

result = train_custom_property_model(
    "my_labels.csv",
    config=CustomPropertyTrainingConfig(
        task_type="regression",
        model="random_forest",
        finetune_method="machine_learning",
        feature_source="sequence",
        mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
        output_dir="runs/random-forest",
        workers=0,
        seed=42,
    ),
)
print(result.metrics["test"])
predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir, trust_model=True,
)
print(predictions.predictions)
```

- **示例结果：** 测试集指标和逐条预测值；具体数值由输入数据与训练过程决定。模型文件保存到 `result.output_dir`。

## 2) 自定义深度学习分类与回归 {#_3}

- **作用：** 用带标签的数据文件训练分类或回归模型，并预测新样本。
- **训练方法：** 11 个登记模型均支持全量微调（`full`）和冻结骨干的线性探测（`linear_probe`）；下例使用 DNABERT-2。模型依赖、表征差异和验证范围见[模型与微调方式](../../faq/tutorial_integrations.md#模型与微调方式)。
- **API：** `dnakit.predictions.train_custom_property_model(source, config=...)`; `dnakit.predictions.predict_custom_properties(source, model_dir=...)`
- **输入：** CSV 的 `sequence`、`target` 列；示例为回归，需先安装 `.[training]`。首次运行可能下载权重；`allow_remote_code=True` 仅用于可信模型。
- **输出：** 返回 `CustomPropertyTrainingResult`，包含训练产物路径和评估指标。
- **示例代码：**

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

result = train_custom_property_model(
    "my_labels.csv",
    config=CustomPropertyTrainingConfig(
        task_type="regression",
        model="dnabert2",
        finetune_method="full",
        allow_remote_code=True,
        mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
        output_dir="runs/my-property",
        workers=0,
        seed=42,
    ),
)
print(result.metrics["test"])
predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir, allow_remote_code=True,
)
print(predictions.predictions)
```

- **示例结果：** 测试集指标和逐条预测值；具体数值由输入数据与训练过程决定。模型文件保存到 `result.output_dir`。

<span id="checkpoint-cache"></span>权重存放位置与自动下载规则见 [FAQ](../../faq/tutorial_integrations.md#property-prediction)。

## 3) NT 任务的全量微调参数

以下配置使用 AdamW、验证 MCC 驱动的 ReduceLROnPlateau，以及连续 5 轮未改善时的早停：

```python
from dnakit.predictions import CustomPropertyTrainingConfig

config = CustomPropertyTrainingConfig(
    task_type="classification",
    model="enformer",
    finetune_method="full",
    learning_rate=1e-4,
    batch_size=64,
    beta1=0.9,
    beta2=0.95,
    weight_decay=0.1,
    gradient_accumulation_steps=1,
    lr_scheduler="reduce_on_plateau",
    scheduler_patience=1,
    scheduler_factor=0.5,
    monitor="mcc",
    monitor_mode="auto",
    early_stopping_patience=5,
    device="cuda:0",
    output_dir="runs/enformer-nt",
)
# 将 config 传给 train_custom_property_model(source, config=config)。
```

`batch_size=64` 是当前选定设备上的每批样本数，`gradient_accumulation_steps=1` 表示每批更新一次。当前训练器在单个设备运行，不自动启动多 GPU 分布式训练。显式配置生效于此次调用；通用默认值仍为 batch size 32、不使用调度器、按验证 loss 早停。

`monitor` 可选分类指标 `mcc`、`f1_macro`、`accuracy`，回归指标 `rmse`、`mae`，或 `validation_loss`。`monitor_mode="auto"` 对分类分数取 `max`，对 loss／回归误差取 `min`，也可显式设置方向。MCC 支持二分类、多分类及多输出分类；多标签分类可用 `f1_macro`。监控的始终是验证集。

每个 epoch 验证结束后调用调度器；按 [PyTorch ReduceLROnPlateau 的 patience 定义](https://docs.pytorch.org/docs/stable/generated/torch.optim.lr_scheduler.ReduceLROnPlateau.html)，`scheduler_patience=1` 容忍一次未改善，连续第二次未改善后将下一轮学习率乘以 0.5。早停的 patience 独立计数，学习率降低不会重置早停计数。最佳 checkpoint、早停和调度器共用监控指标及绝对改善阈值 `monitor_min_delta`（默认 `1e-12`，调度器使用 `threshold_mode="abs"`），可正确处理零或负 MCC。训练结束恢复最佳监控分数对应的权重。

`history.json` 保存 `monitor_value`、`primary_metric_name`、本轮实际使用的 `learning_rate` 和调度后的 `next_learning_rate`；`model.json` 的 `training` 保存完整配置。`epochs` 决定训练轮数上限。相同参数也可通过 `dnakit train-property` 的 `--monitor`、`--lr-scheduler`、`--scheduler-patience`、`--scheduler-factor`、`--early-stopping-patience` 等选项传入。

## 4) 折训练与完整训练集重训

神经网络的 `full`／`linear_probe` 训练允许省略 validation 或 test 分区。使用 `DataMapping(split_field="split")` 时，可以只提供 `train`／`valid` 两种行；训练器只计算存在的分区，不会为折训练额外拆出 test。有验证集时输出 `validation_predictions.jsonl`，其中包含 `record_id`、预测值和各类别概率；`test_predictions.jsonl` 仅在实际提供 test 时生成。

设置 `log_train_metrics=True` 可在每轮结束后额外评估训练集，用于绘制 train／valid MCC 曲线；这会增加一次训练集前向计算。`history.json` 包含 `train_mcc`、`validation_mcc`，不适用或未计算的值为 `null`。`train_loss` 是该轮训练批次损失按样本数加权的平均值，验证 loss 和两侧 MCC 是该轮结束时评估的结果。

选定超参数与轮数后，可以仅使用完整 train 重训：

```python
from dnakit.predictions import CustomPropertyTrainingConfig, train_custom_property_model

config = CustomPropertyTrainingConfig(
    task_type="classification", model="enformer", finetune_method="full",
    checkpoint_path="models/enformer", output_dir="runs/enformer-final",
    train_ratio=1.0, validation_ratio=0.0, test_ratio=0.0,
    epochs=8,  # 事先由训练集内部验证确定的轮数
    learning_rate=1e-4, batch_size=64, beta1=0.9, beta2=0.95, weight_decay=0.1,
    gradient_accumulation_steps=1,
    early_stopping_patience=0, lr_scheduler="none", log_train_metrics=True,
)
result = train_custom_property_model("official_train.csv", config=config)
```

没有验证集时必须关闭早停和验证驱动调度；训练执行固定轮数，保存最后一轮权重，日志以 `train_loss` 记录进展。官方 test 保持独立，使用最终模型通过 `predict_custom_properties` 单独推理。`model.json` 保存实际划分与配置，返回指标不包含不存在的分区。机器学习模式未启用 CV 时仍要求非空 train／validation／test；神经网络多折训练需在用户脚本中显式循环，不能使用仅适用于机器学习模式的 `cv_folds`。
