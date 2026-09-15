# 性质预测：深度学习与机器学习

新增的 AlphaGenome、ESM、Evo 2、NVIDIA NIM 等 **14 个已有模型调用接口**见[已训练模型的预测接口](../native/models.md)，每个函数均列出参数、返回值及已登记的凭据要求。CADD、AlphaMissense、SpliceAI 等变异注释与评分入口见[变异分析接口](../native/variants.md)。下方继续介绍 DNAKit 的自定义训练与预测流程。

缺少或无法导入预测依赖时，DNAKit 会先向终端标准错误输出安装提示，再抛出原有异常；异常中也保留安装命令。`verbose=False` 关闭运行进度，依赖错误提示仍保留。DNAKit 不会自动执行安装。

在 **DNAKit 项目根目录**，按所用功能安装：

| 功能 | 安装命令 |
| --- | --- |
| Enformer（含 27 个分类任务） | `python -m pip install -e ".[neural,neural-enformer]"` |
| 通用 PyTorch / Transformers 主干 | `python -m pip install -e ".[neural]"` |
| 自定义深度学习训练 | `python -m pip install -e ".[training]"` |
| 基础机器学习 | `python -m pip install -e ".[ml]"` |
| XGBoost / LightGBM / CatBoost | 分别使用 `.[xgboost]` / `.[lightgbm]` / `.[catboost]` |

通过发布包安装的用户，可将 `-e ".[选项]"` 改成 `"dnakit[选项]"`，无需源码目录。Evo 2、Caduceus 及需要官方源码的模型会显示各自的安装提示；Evo 2 的 extra 要求 Python 3.11–3.12，AlphaGenome、LucaOneTasks、JanusDNA 仍需满足官方环境要求。

多个模型可以共用满足各自版本约束的 PyTorch 或 TensorFlow；二者不能互相替代。上述命令不强制升级或重新安装已满足要求的包，只补齐缺失或不兼容的依赖。当前 `neural` / `training` 的 PyTorch 范围是 `>=2.5,<3`，没有为每个模型固定不同版本；第三方模型仍可能有更严格的约束。参见 [pip 安装行为](https://pip.pypa.io/en/stable/cli/pip_install/)。

本页包含三条路径：

1. **用户自定义性质：** 从用户文件训练分类或回归模型，支持全量微调和线性探测；
2. **官方现成性质：** 加载官方已经发布的任务头直接推理，不在本地重新训练。
3. **机器学习性质：** 使用 scikit-learn 或可选的 XGBoost、LightGBM、CatBoost，训练自己的分类或回归模型。

原有的[神经网络表征](08_fingerprints.md)仍只负责提取 embedding。下面先介绍自定义训练；随后编号 1–54 的内容仍表示 54 个官方直接推理功能。

本页当前共 54 个可直接推理的功能。下面每个功能都单独说明作用、预测方法、API、输入与输出。示例结果中的占位符只表示输出结构；具体数值、标签数和轨道数取决于输入、物种、组织筛选条件及实际 checkpoint。

本章涉及的模型论文、官方仓库和权重来源统一收录在 [FAQ：深度学习性质预测的参考文献](../../faq.md#deep-learning-property-prediction-references)。

## 输入与输出的异常处理

直接预测拒绝空序列和非法字符，并在序列错误中提供 `record_id`。Enformer 任务允许
1、199、200、201 bp 等长度；默认上限为 196608 bp，显式 `max_length` 可以收紧限制，
不能放宽模型上限。不会静默截断输入；IUPAC 符号按 `ambiguity_policy` 处理。

`PredictionOutput` 拒绝空、非有限、复数或无法转换的数值，以及重复输出列名。
直接分类／分割任务要求输出为 `[0,1]` 概率；Enformer 分类还要求输出维度等于类别数、
概率和在 `1e-5` 容差内等于 1。连续分数、回归和轨道输出不套用概率范围。
自定义后端也必须遵守所选任务的输出约定。

直接预测及 embedding 公共入口将内存不足（包括后端异常链中的 OOM）报告为
`MODEL_OUT_OF_MEMORY`，提示减小批次、释放内存或使用支持的其他设备；不自动输出替代预测。
显式指定不存在的 CUDA 设备或编号会在模型加载前报错。

## 权重缓存与自动下载 {#checkpoint-cache}

预测前会检查本地预训练权重：缺失时只下载当前模型或任务所需的文件，并提示开始、完成及本地绝对路径；已有权重直接复用。`verbose=False` 可关闭提示。Enformer 任务默认存放在 `./ckpt/enformer-benchmarks/`，其他主干存放在 `./ckpt/<模型目录>/`。机器学习使用的 embedding 主干同样按需下载。用户训练的分类器（`estimator.joblib`）和神经网络任务头（`head.npz`）没有登记公共下载来源，需保留 `model_dir` 中的训练输出。LucaOneTasks 使用所提供源码目录中的官方下载流程。

Enformer 的 27 个分类任务无需指定权重路径即可预测。使用 `checkpoint_dir` 可指定缓存目录：已有文件也会在其 `enformer-benchmarks/` 及对应任务类别的 `nt-revised/` 或 `genomic-benchmarks/` 子目录中查找；缺失时直接下载到指定目录。也可使用 `checkpoint_path` 指定精确文件路径，所选 Enformer 任务的文件缺失时会下载到该路径。两个参数不能同时使用。

安装命令为 `python -m pip install "dnakit[neural,neural-enformer]"`，其中已包含 Google Drive 下载所需的 `gdown`。下载先写入临时文件，检查后再放入正式路径；失败或中断的文件不会被当作完整权重。预测时继续校验任务和模型参数，类别含义以 `output.output_names` 中的标签为准。

## 机器学习分类与回归

与深度学习共用 `CustomPropertyTrainingConfig`、`train_custom_property_model()` 和 `predict_custom_properties()`。设置 `finetune_method="machine_learning"` 后，`model` 表示机器学习算法；embedding 主干另用 `embedding_model` 指定。

安装基础机器学习及表格读取依赖（不需要 PyTorch）：

```bash
python -m pip install -e ".[ml]"
```

安装三个可选提升树库：

```bash
python -m pip install -e ".[xgboost,lightgbm,catboost]"
```

如需 DNA 模型 embedding，再安装 `.[neural]` 及所选主干要求的附加依赖。缺少某个可选算法库时只会影响该算法，调用时会给出对应安装命令。

### 模型覆盖

目前登记 37 类算法。以下是完整模型名称；`available_machine_learning_models()` 返回每类支持的任务和依赖，不会导入或安装模型库。

| 任务范围 | `model` |
| --- | --- |
| 分类和回归 | `ridge`、`sgd`、`svm`、`linear_svm`、`nu_svm`、`knn`、`decision_tree`、`extra_tree`、`random_forest`、`extra_trees`、`gradient_boosting`、`hist_gradient_boosting`、`adaboost`、`bagging`、`gaussian_process`、`mlp`、`dummy`、`xgboost`、`lightgbm`、`catboost` |
| 仅分类 | `logistic_regression`、`perceptron`、`gaussian_nb`、`multinomial_nb`、`bernoulli_nb`、`complement_nb`、`lda`、`qda` |
| 仅回归 | `linear_regression`、`lasso`、`elastic_net`、`bayesian_ridge`、`ard`、`huber`、`ransac`、`theil_sen`、`kernel_ridge` |

`mlp` 是 scikit-learn 的浅层神经网络，`dummy` 是基线。这里的范围是监督性质预测，不包括聚类、异常检测及自动模型搜索。各算法仍有自己的样本量、参数和计算成本约束；例如高斯过程不适合直接用于超大样本集。

### 三种输入特征

| `feature_source` | 行为 |
| --- | --- |
| `sequence`（默认） | 序列长度、GC 比例、模糊碱基比例，加固定顺序的 k-mer 频率；`kmer_k=3`，允许 1–6，维度为 `3 + 4**kmer_k`。含模糊碱基的窗口跳过，GC 分母只含 A/C/G/T；无有效窗口或碱基时相应频率为 0。无需深度学习权重。 |
| `embedding` | 冻结 `embedding_model`（默认 `dnabert2`）提取向量，再训练机器学习模型。沿用 checkpoint、pooling、batch_size、device 等表征参数。 |
| `numeric` | 使用显式 `DataMapping.feature_fields`，按给出的列顺序读取数值；可设置 `sequence_field=None`，输入文件不需要序列。缺失、无穷、非数值特征会报错，当前不自动填补。 |

数值特征列必须显式指定，不能与标签、ID、split 或序列字段重叠。推理时必须保留训练时的特征列名称与顺序。不要把目标派生列或其他泄漏信息选为特征。

下面使用自己的序列数据训练随机森林回归：

```python
from dnakit.predictions import (
    CustomPropertyTrainingConfig, DataMapping,
    train_custom_property_model, predict_custom_properties,
)

config = CustomPropertyTrainingConfig(
    task_type="regression",  # 分类改为 classification
    finetune_method="machine_learning",
    model="random_forest",
    model_params={"n_estimators": 300, "max_depth": 12},
    feature_source="sequence",
    kmer_k=3,
    mapping=DataMapping(sequence_field="sequence", label_fields=("target",)),
    output_dir="runs/random-forest",
    workers=6,
    seed=42,
)
result = train_custom_property_model("my_labels.csv", config=config)
print(result.metrics["test"])

predictions = predict_custom_properties(
    "new_sequences.csv", model_dir=result.output_dir,
    trust_model=True,  # 仅加载自己训练或来源可信的模型
)
print(predictions.predictions)
```

使用数值特征时改为 `feature_source="numeric"`，并传入例如 `DataMapping(sequence_field=None, feature_fields=("gc", "length"), label_fields=("target",))`；使用 embedding 时改为 `feature_source="embedding", embedding_model="dnabert2"`，并按主干要求设置 `allow_remote_code` 等参数。

### 训练、评估与保存

默认仍为种子 42、80%/10%/10% 划分；已映射的 split 列优先。标签模式、测试指标和深度学习路径一致。多标签、多输出分类及多目标回归按目标分别训练估计器；训练中恒定的二值标签使用恒定分类器。分类预测取最大概率类别，多标签使用 `threshold`。

`model_params` 直接传给所选估计器构造函数，接受可序列化的参数，不接受任意 Python 估计器对象。`workers` 控制默认 `n_jobs` 和数值库线程数（0 按 1 处理），CatBoost 使用 `thread_count`；模型参数可覆盖估计器的并行数和随机种子。`standardize_features=True` 默认在训练集上拟合标准化器，验证和测试只做变换。`multinomial_nb`、`complement_nb` 要求非负输入，使用保持非负性的 MaxAbsScaler；不会把负 embedding 静默平移成非负数据。

无 `predict_proba` 的分类器（如 RidgeClassifier、LinearSVC）使用训练集内部最多 3 折 sigmoid 概率校准，每类至少需要 2 个训练样本。标准化也在每个校准训练折内拟合，验证集和测试集不参与。校准行为记录在结果 warnings 中。参见 [scikit-learn 概率校准 API](https://scikit-learn.org/stable/modules/generated/sklearn.calibration.CalibratedClassifierCV.html)。

机器学习路径不使用深度学习的 `epochs`、AdamW 参数或梯度累积；迭代次数、树数量等通过 `model_params` 设置。验证集只用于评估，不会自动调参或选择最优模型。

回归标签沿用现有接口的训练集均值/标准差归一化，报告和预测恢复到原始单位。分类 `class_weight` / `class_weights` 字典使用从 0 开始的编码类别索引；类别顺序保存在 `model.json` 的 `target_schema.class_names` 中，也可使用算法支持的 `"balanced"` 等选项。

输出包括 `model.json`、`estimator.joblib`、`metrics.json`、`data_inspection.json`、`splits.csv`、`test_predictions.jsonl`。估计器文件包含拟合好的特征变换；manifest 记录特征顺序、配置、标签语义和 Python/依赖版本。机器学习没有神经网络训练 history，结果的 `history` 为空。

joblib 加载可能执行 Python 代码，因此机器学习推理必须显式设置 `trust_model=True`；仅用于自己训练或来源可信的工件。重新加载会核对 manifest 记录的 Python 和依赖版本，版本不一致时要求使用原训练环境。该约束依据 [scikit-learn 模型持久化说明](https://scikit-learn.org/stable/model_persistence.html)。

可直接编辑并运行 `others/examples/train_machine_learning_properties.py`；`--demo` 使用人工数据验证完整训练和重新加载流程：

```bash
python others/examples/train_machine_learning_properties.py --demo
```

## 用户自定义深度学习分类与回归

先安装训练与常见表格格式依赖：`python -m pip install "dnakit[training]"`。

`train_custom_property_model()` 接受一个本地数据文件，自动或按用户指定的映射找出 ID、DNA/RNA 序列、一个或多个标签以及可选 split。默认随机种子为 42，训练/验证/测试比例为 0.8/0.1/0.1；若映射了已有 split 列，则使用文件中的划分。`task_type` 必须由用户明确指定为 `classification` 或 `regression`，不会仅凭标签外观替用户决定科学任务。

内置读取器支持 CSV、TSV、其他单字符分隔文本、JSON、JSONL、YAML、XML、Excel、Parquet、SQLite、FASTA、FASTQ 和 GenBank。嵌套 JSON/YAML/XML 可使用点路径；FASTA/FASTQ 描述中的 `key=value` 会映射为 `header.key`。其他格式可实现 `DataReader` 协议。大模型只生成受约束的声明式 `DataMapping`，不会生成或执行 Python、SQL、正则表达式或 shell 代码。

### 标签语义

- 一个分类标签列：二分类或多分类；
- 多个二值标签列：自动识别 one-hot 多分类或多标签分类，也可用 `classification_mode` 明确指定；
- 多个非二值分类列：多输出分类，每列保留独立类别表；
- 一个或多个数值标签列：单目标或多目标回归。

分类报告包括 Accuracy、Balanced Accuracy、Precision、Recall、F1、Jaccard、MCC、Cohen's kappa、Log loss、AUROC、Average Precision、混淆矩阵及逐类别/逐标签结果；多标签任务另含 subset accuracy、Hamming loss、LRAP 和 coverage error。回归报告包括 MAE、MSE、RMSE、R²、explained variance、median absolute error、max error、bias、MAPE、SMAPE、Pearson 和 Spearman，并同时给出汇总与逐目标结果。F1、AUROC 等分类指标不会伪装成回归指标；样本不足或数学上未定义的指标返回 `null`。

### 模型与微调方式

默认 `finetune_method="full"`。当前 Transformers 路径支持 DNABERT-2、NT-v2、HyenaDNA、Caduceus、GROVER、LucaOne 和 GENERator 的全量微调；全部已登记表征模型均可做 `linear_probe`，包括 AlphaGenome、Enformer、Evo 2 和 JanusDNA。请用 `available_custom_training_models()` 获取当前版本的实际支持矩阵。需要 checkpoint 自带代码的模型仍必须显式设置 `allow_remote_code=True`。

### 训练参数（Python API）

直接在 `CustomPropertyTrainingConfig` 中修改参数，无需命令行参数。全量微调与线性探测使用相同默认值：

| Python 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `workers` | `6` | CPU 数据加载进程数；`0` 表示在主进程加载 |
| `batch_size` | `32` | 每个训练设备的单批样本数，不是梯度累积后的批大小 |
| `learning_rate` | `1e-4` | AdamW 学习率，即 lr |
| `beta1` | `0.9` | AdamW 一阶矩衰减系数 |
| `beta2` | `0.95` | AdamW 二阶矩衰减系数 |
| `weight_decay` | `0.1` | AdamW 权重衰减 |

`workers` 用于训练、验证和测试的数据取样及 CPU 批组装，编码器和 GPU 运算仍在主进程；它不控制初次文件解析或大模型请求的并发数。固定种子决定划分和每轮样本顺序，调整 workers 不会改变顺序。当前训练器在 `device` 指定的单个设备运行（例如 `"cuda:0"`），尚未实现多卡 DDP；不会根据显卡数量自动增大批大小。`gradient_accumulation_steps` 默认 `1`，可通过 Python 配置。

使用 `workers>0` 的 Python 脚本需要下例的 `if __name__ == "__main__":` 保护；Notebook 或不支持子进程的环境建议设置 `workers=0`。数据加载进程在训练完成、提前停止或异常时释放。

```python
import os

from dnakit.predictions import (
    CustomPropertyTrainingConfig,
    LLMSchemaConfig,
    train_custom_property_model,
)

def main():
    result = train_custom_property_model(
        "my_labels.xlsx",
        config=CustomPropertyTrainingConfig(
            task_type="regression",  # 或 "classification"
            model="dnabert2",
            finetune_method="full",  # 或 "linear_probe"
            output_dir="runs/my-property",
            workers=6,
            batch_size=32,
            learning_rate=1e-4,
            beta1=0.9,
            beta2=0.95,
            weight_decay=0.1,
            seed=42,
            train_ratio=0.8,
            validation_ratio=0.1,
            test_ratio=0.1,
            allow_remote_code=True,
            llm=LLMSchemaConfig(
                provider="deepseek",
                model="deepseek-v4-flash",
                api_key=os.environ["DEEPSEEK_API_KEY"],
            ),  # 不需要大模型辅助时省略 llm
        ),
    )
    print(result.metrics["test"])


if __name__ == "__main__":
    main()
```

训练结果目录包含 `model.json`、`head.npz`、`metrics.json`、`history.json`、`data_inspection.json`、`splits.csv` 和 `test_predictions.jsonl`；全量微调还包含 `backbone/`。`model.json` 的 `training` 字段记录实际 workers、批大小、学习率、两个 beta 和 weight decay。`predict_custom_properties()` 可重新加载这些产物进行推理。

仓库提供可直接编辑参数的脚本 `others/examples/train_custom_properties.py`：修改输入文件与配置后，在 IDE 中运行即可；脚本默认不启用大模型辅助。

### 大模型服务商与数据分享范围

不仅支持 DeepSeek。`LLMSchemaConfig(provider=..., model=...)` 内置以下连接预设，统一调用 OpenAI 兼容的 `POST /chat/completions`，不需要安装各厂商 SDK。`available_schema_llm_providers()` 可离线查询这些预设。除兼容旧用法的 DeepSeek 默认模型外，必须自行填写账号可用的模型 ID；配置不会自动切换模型或服务商。

| `provider` | 默认 `base_url` | 默认密钥环境变量 |
| --- | --- | --- |
| `deepseek` | `https://api.deepseek.com` | `DEEPSEEK_API_KEY` |
| `gemini` | `https://generativelanguage.googleapis.com/v1beta/openai` | `GEMINI_API_KEY` |
| `qwen` | `https://dashscope.aliyuncs.com/compatible-mode/v1` | `DASHSCOPE_API_KEY` |
| `kimi` | `https://api.moonshot.cn/v1` | `MOONSHOT_API_KEY` |
| `openrouter` | `https://openrouter.ai/api/v1` | `OPENROUTER_API_KEY` |
| `openai` | `https://api.openai.com/v1` | `OPENAI_API_KEY` |
| `custom` | 必须显式指定 | `LLM_API_KEY` |

例如，仅分享字段与类型结构，使用 Google Gemini 辅助分析：

```python
from dnakit.predictions import LLMSchemaConfig, inspect_property_data

llm = LLMSchemaConfig(
    provider="gemini",
    model="gemini-3.8-flash",  # 可替换为账号可用的文本/JSON 模型。
    include_sample_values=False,
)
inspection = inspect_property_data("my_labels.xlsx", llm=llm)
print(inspection.mapping.to_dict())
print(inspection.recommendations)
# 训练时将同一个 llm 传给 CustomPropertyTrainingConfig(llm=llm, ...)。
```

`api_key` 仍可在 Python 中传入，或用 `api_key_env` 指定自己的环境变量名；不要将真实密钥写进示例或版本库。`base_url` 可覆盖预设，密钥只发送到你配置的地址，HTTP 重定向被拒绝。Gemini 此处使用其[官方 OpenAI 兼容接口](https://ai.google.dev/gemini-api/docs/openai)，不是原生 `generateContent` 或 Vertex AI 接口。

千问默认地址为北京地域，官方仍支持该域名；也可以填控制台提供的业务空间专属地址，例如 `https://<业务空间ID>.cn-beijing.maas.aliyuncs.com/compatible-mode/v1`。地域地址和 API Key 必须匹配，见[千问接口说明](https://help.aliyun.com/zh/model-studio/compatibility-of-openai-with-dashscope)。Kimi 国际站可使用 `base_url="https://api.moonshot.ai/v1"`，参见 [Kimi 快速开始](https://platform.kimi.ai/docs/overview)。OpenRouter 的模型名须使用其平台的完整 ID，且应检查目标模型支持的参数，见 [OpenRouter 快速开始](https://openrouter.ai/docs/quickstart)和[结构化输出说明](https://openrouter.ai/docs/guides/features/structured-outputs)；数据将由其平台路由到所选模型的提供方。

兼容参数可以直接在 Python 中调整：

- `temperature=None`（默认）：不发送温度参数，适用于不接受自定义温度的模型；需要时可显式设置。
- `token_limit_parameter="auto"`（默认）：OpenAI 使用 `max_completion_tokens`，其他兼容服务使用 `max_tokens`；也可显式选择这两个字段之一。依据 [OpenAI Docs 的 Chat Completions 参数说明](https://developers.openai.com/api/reference/python/resources/chat/subresources/completions/methods/create)，旧 `max_tokens` 不兼容部分推理模型。
- `json_mode=True`（默认）：要求 JSON Object 输出。若模型不支持该接口参数，可设为 `False`；提示词仍要求 JSON，返回结果仍须经过严格映射校验。
- `extra_body`：传递模型专用参数，例如支持关闭思考的千问模型可设 `{"enable_thinking": False}`；支持该选项的 Kimi 模型可设 `{"thinking": {"type": "disabled"}}`。不会统一强行关闭所有模型的思考模式，需按[千问](https://help.aliyun.com/zh/model-studio/qwen-structured-output)或 [Kimi 模型参数](https://platform.kimi.ai/docs/api/models-overview)确认支持情况。请求当前为非流式；只支持流式思考的模型需关闭思考或改用支持非流式的模型。

`extra_body` 不允许覆盖消息、鉴权信息、模型、输出 token 上限等核心字段。输出被截断时会明确报错，应调整 `max_output_tokens` 或模型专用思考参数，不会自动扩大额度重试。服务商预设表示接口接入能力，不保证该平台的每个模型都支持此任务；当前已做模拟接口测试，未用各平台付费密钥逐一联网验证。

默认最多分享 8 条样本、12,000 个字符及结构摘要，**不会整表上传**；小文件也至少保留一条完整记录不分享。`include_sample_values=False` 只发送字段/类型结构。`DataInspection` 记录分享的样本数和字符数；密钥不进入请求正文、结果、异常展示或训练产物。

## 1) RNA-seq 覆盖与表达

- **作用：** 根据 DNA 序列预测 RNA-seq 覆盖和基因表达轨道。
- **预测方法：** 使用 AlphaGenome all-folds 官方权重直接推理 `RNA_SEQ` 输出；短序列会居中补 `N` 到官方支持的上下文长度，可用 `ontology_terms` 限定组织或细胞类型。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="rna_seq")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也可改用 `VariantContext` 和 `predict_variant_effects()`，得到 REF、ALT、ALT−REF 三组轨道。支持人或小鼠。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="rna_seq",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <RNA-seq 轨道数>)
('position', 'track')
```

## 2) CAGE 转录起始与表达

- **作用：** 预测 CAGE 信号，用于描述加帽转录本的转录起始位置和表达强度。
- **预测方法：** 使用 AlphaGenome 官方 `CAGE` 输出头直接生成位置 × 轨道矩阵；`ontology_terms` 可筛选目标组织或细胞类型。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="cage")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext` 形式的单 SNV 上下文。序列最长 1,048,576 bp，支持人或小鼠。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="cage",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <CAGE 轨道数>)
('position', 'track')
```

## 3) PRO-cap 新生转录起始

- **作用：** 预测 PRO-cap 新生 RNA 5′ 端信号，用于定位活跃转录起始事件。
- **预测方法：** 使用 AlphaGenome 官方 `PROCAP` 输出头直接推理；输出分辨率记录在 `metadata["resolution_bp"]` 中。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="procap")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可指定 `organism` 和 `ontology_terms`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="procap",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["resolution_bp"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <PRO-cap 轨道数>)
<每个输出位置对应的 bp 数>
```

## 4) ATAC-seq 染色质可及性

- **作用：** 预测 ATAC-seq 开放染色质轨道，用于识别可能处于可及状态的调控区域。
- **预测方法：** 使用 AlphaGenome 官方 `ATAC` 输出头直接推理，并按位置返回不同组织或细胞类型的可及性信号。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="atac")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持只包含一个碱基变化的 `VariantContext`。可通过 `ontology_terms` 限制轨道。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="atac",
        ontology_terms=("UBERON:0002048",),
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names[:2])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <筛选后的 ATAC-seq 轨道数>)
(<前两个官方轨道名称>,)
```

## 5) DNase-seq 染色质可及性

- **作用：** 预测 DNase-seq 超敏信号，用于描述染色质可及性。
- **预测方法：** 使用 AlphaGenome 官方 `DNASE` 输出头直接推理；短输入会居中补齐，补齐长度保存在结果元数据中。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="dnase")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可选 `organism="human"` 或 `"mouse"`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="dnase",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["input_length"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <DNase-seq 轨道数>)
16384
```

## 6) 组蛋白修饰

- **作用：** 预测组蛋白 ChIP-seq 轨道，用于表征不同组蛋白修饰在序列区域中的信号。
- **预测方法：** 使用 AlphaGenome 官方 `CHIP_HISTONE` 输出头直接推理，返回位置 × 组蛋白轨道矩阵。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="chip_histone")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可按物种和 ontology term 选择输出。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="chip_histone",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <组蛋白 ChIP-seq 轨道数>)
('position', 'track')
```

## 7) 转录因子结合

- **作用：** 预测转录因子 ChIP-seq 结合轨道，用于定位可能的转录因子结合信号。
- **预测方法：** 使用 AlphaGenome 官方 `CHIP_TF` 输出头直接推理，返回位置 × 转录因子轨道矩阵。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="chip_tf")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可用 `ontology_terms` 限定组织或细胞类型。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="chip_tf",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names[:2])
```

- **示例结果（结构示意）：**

```text
(<位置或 bin 数>, <TF ChIP-seq 轨道数>)
(<前两个官方轨道名称>,)
```

## 8) 剪接供体/受体位点

- **作用：** 预测序列各位置的剪接供体和剪接受体相关轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_SITES` 输出头直接推理，并保留输出分辨率及轨道名称。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_sites")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`，变异结果额外包含 ALT−REF 差值轴。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="splice_sites",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<位置数>, <供体/受体轨道数>)
('position', 'track')
```

## 9) 剪接位点使用比例

- **作用：** 预测候选剪接位点的使用信号或使用比例轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_SITE_USAGE` 输出头直接推理；具体轨道身份由 `output_names` 给出。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_site_usage")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`。可指定物种、ontology term 和官方上下文长度。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="splice_site_usage",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["resolution_bp"])
```

- **示例结果（结构示意）：**

```text
(<位置数>, <剪接位点使用轨道数>)
<每个输出位置对应的 bp 数>
```

## 10) 剪接连接计数

- **作用：** 预测候选 splice junction 的连接计数轨道。
- **预测方法：** 使用 AlphaGenome 官方 `SPLICE_JUNCTIONS` 输出头直接推理；结果元数据中的 `junctions` 保存候选连接标识。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="splice_junctions")[必须], backend[可选])`
- **输入：** 只接受 `BiologicalSequence`。该任务的候选连接集合可能随等位基因变化，因此当前不接受 `VariantContext`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="splice_junctions",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
print(len(output.metadata.get("junctions", ())))
```

- **示例结果（结构示意）：**

```text
(<候选剪接连接数>, <轨道数>)
('junction', 'track')
<候选剪接连接数>
```

## 11) 三维 DNA 接触图

- **作用：** 预测染色质区域之间的三维 DNA 接触概率图。
- **预测方法：** 使用 AlphaGenome 官方 `CONTACT_MAPS` 输出头直接推理，输出两个位置轴和一个轨道轴。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="alphagenome", task="contact_maps")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持 `VariantContext`，变异结果在最前面增加 REF、ALT、ALT−REF 轴。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="alphagenome",
        task="contact_maps",
        organism="human",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["axes"])
```

- **示例结果（结构示意）：**

```text
(<bin 数>, <bin 数>, <接触图轨道数>)
('position_1', 'position_2', 'track')
```

## 12) 人类调控轨道

- **作用：** 预测人类基因组的 5,313 条 CAGE、DNase/ATAC、转录因子 ChIP 和组蛋白 ChIP 轨道。
- **预测方法：** 使用公开的 Enformer 预训练权重直接推理；默认把输入居中补 `N` 到 196,608 bp，并输出 896 个 128 bp bin。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="enformer", task="human_tracks")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也可用 `VariantContext` 和 `predict_variant_effects()`，获得 REF、ALT、ALT−REF。`max_length` 必须是 114,688–393,216 范围内的 128 倍数。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="enformer",
        task="human_tracks",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["bin_size_bp"])
```

- **示例结果（默认上下文）：**

```text
(896, 5313)
128
```

## 13) 小鼠调控轨道

- **作用：** 预测小鼠基因组的 1,643 条 CAGE、DNase/ATAC、转录因子 ChIP 和组蛋白 ChIP 轨道。
- **预测方法：** 使用 Enformer 的小鼠输出头直接推理；默认上下文为 196,608 bp，输出为 128 bp 分箱轨道。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="enformer", task="mouse_tracks")[必须], backend[可选])`
- **输入：** `BiologicalSequence`；也支持只含一个 SNV 的 `VariantContext`。输入不足上下文长度时自动在两端居中补 `N`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("mouse-region-1", "ACGT" * 4096)],
    config=PropertyPredictionConfig(
        model="enformer",
        task="mouse_tracks",
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["organism"])
```

- **示例结果（默认上下文）：**

```text
(896, 1643)
mouse
```

## 14) 14 类单碱基基因组分割

- **作用：** 为人 DNA 序列中的每个碱基预测 14 类基因和调控元件概率。
- **预测方法：** 使用带训练完成分割头的 `InstaDeepAI/segment_nt`，即 NT-v2 骨干网络加 U-Net 分割头；对二分类 logits 做 softmax 后返回“该特征存在”的概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="segmentnt", task="genomic_segmentation", allow_remote_code=True)[必须], backend[可选])`
- **输入：** 最长 30,000 bp 的人 `BiologicalSequence`。14 个类别分别列在下方的 14.1–14.14 子节中。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("region-1", "ACGT" * 1000)],
    config=PropertyPredictionConfig(
        model="segmentnt",
        task="genomic_segmentation",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(len(output.output_names))
```

- **示例结果：**

```text
(4000, 14)
14
```

以下 14 个类别共用上述预测接口，每类对应 `output.values` 的一列，实际列名与顺序以 `output.output_names` 为准。

### 14.1) 蛋白质编码基因（protein-coding gene）

预测每个碱基属于蛋白质编码基因区域的概率。

### 14.2) 长链非编码 RNA（lncRNA）

预测每个碱基属于长链非编码 RNA 基因区域的概率。

### 14.3) 外显子（exon）

预测每个碱基属于外显子区域的概率。

### 14.4) 内含子（intron）

预测每个碱基属于内含子区域的概率。

### 14.5) 剪接供体位点（splice donor）

预测每个碱基属于剪接供体位点的概率。

### 14.6) 剪接受体位点（splice acceptor）

预测每个碱基属于剪接受体位点的概率。

### 14.7) 5′ 非翻译区（5′ UTR）

预测每个碱基属于 5′ 非翻译区的概率。

### 14.8) 3′ 非翻译区（3′ UTR）

预测每个碱基属于 3′ 非翻译区的概率。

### 14.9) CTCF 结合区域（CTCF-bound）

预测每个碱基属于 CTCF 结合区域的概率。

### 14.10) 多聚腺苷酸化信号（polyA signal）

预测每个碱基属于多聚腺苷酸化信号的概率。

### 14.11) 组织特异性增强子（tissue-specific enhancer）

预测每个碱基属于组织特异性增强子区域的概率。

### 14.12) 组织不变性增强子（tissue-invariant enhancer）

预测每个碱基属于组织不变性增强子区域的概率。

### 14.13) 组织特异性启动子（tissue-specific promoter）

预测每个碱基属于组织特异性启动子区域的概率。

### 14.14) 组织不变性启动子（tissue-invariant promoter）

预测每个碱基属于组织不变性启动子区域的概率。

## 15) 长上下文零样本变异效应

- **作用：** 计算单 SNV 的 REF、ALT 序列似然及 ALT−REF 变异效应分数。
- **预测方法：** 使用 Evo 2 7B 官方基础模型分别计算完整 REF 和 ALT 上下文的平均序列 log-likelihood，再计算 `alternate − reference`。该值是零样本排序分数，不是临床校准的致病概率。
- **API：** `dnakit.predictions.predict_variant_effects(inputs[必须], config=PropertyPredictionConfig(model="evo2", task="variant_effect")[必须], backend[可选])`
- **输入：** `VariantContext`；REF 与 ALT 必须等长且只有一个 A/C/G/T 位点不同。
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

## 16) 外显子概率

- **作用：** 根据一个候选位点两侧的正向和反向上下文预测该位点属于 exon 的概率。
- **预测方法：** 从 Evo 2 `blocks.26` 读取正向、反向上下文的末 token 4,096 维 embedding，拼接为 8,192 维，再通过官方已经训练好的 MLP 头输出概率；不会重新训练分类器。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="evo2", task="exon_probability", allow_remote_code=True)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`，顺序必须是 forward context、reverse context；两条都使用 `sequence_type="gene"`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

contexts = BiologicalSequencePair(
    "site-1",
    BiologicalSequence("forward-1", "ACGT" * 256),
    BiologicalSequence("reverse-1", "ACGT" * 256),
)
result = predict_pair_properties(
    [contexts],
    config=PropertyPredictionConfig(
        model="evo2",
        task="exon_probability",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签取决于概率和阈值）：**

```text
('exon_probability',)
exon 或 non_exon
```

## 17) 等位基因条件概率变异效应

- **作用：** 根据变异位点上游序列，比较下一 token 支持 REF 与 ALT 碱基的概率。
- **预测方法：** 使用 GENERator v2 eukaryote 1.2B 官方权重；最多读取 8,192 bp 上游上下文，去除开头 `N` 并左截断到 6 的倍数，然后计算 `log(p_ref/(p_alt+1e-10))`。该值是零样本排序分数，不是临床校准概率。
- **API：** `dnakit.predictions.predict_variant_effects(inputs[必须], config=PropertyPredictionConfig(model="generator", task="variant_effect", allow_remote_code=True)[必须], backend[可选])`
- **输入：** `VariantContext`；REF 与 ALT 必须等长且只含一个 SNV。模型只使用变异位点之前的参考上下文。
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

## 18) 核酸–蛋白中心法则关系

- **作用：** 预测一条 gene 序列与一条 protein 序列是否具有 CentralDogma 关系。
- **预测方法：** 调用 LucaOneTasks 官方 `CentralDogma` 二分类 checkpoint 和 `src/predict_v1.py` 直接推理，返回正类概率及官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="central_dogma", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；第一条必须为 `sequence_type="gene"`，第二条必须为 `sequence_type="protein"`。`model_source_path` 必须指向 LucaOneTasks 官方源码目录。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "central-dogma-1",
    BiologicalSequence("gene-1", "ATGGCC" * 100, "gene"),
    BiologicalSequence("protein-1", "MAMAPRTEINSTRING" * 10, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="central_dogma",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```

## 19) SupKTax 分类

- **作用：** 使用官方 SupKTax checkpoint 对 gene 序列进行分类。
- **预测方法：** 调用 LucaOneTasks 官方 `SupKTax` 多分类 checkpoint 直接推理，并返回概率最高的 `top_k` 个 `label.txt` 标签。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="supktax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`。结果标签以该 checkpoint 的 `label.txt` 为准。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="supktax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个官方标签>,)
```

## 20) GenusTax 分类

- **作用：** 使用官方 GenusTax checkpoint 对 gene 序列进行分类。
- **预测方法：** 调用 LucaOneTasks 官方 `GenusTax` 多分类 checkpoint 直接推理，并按概率返回 Top-k 标签。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="genustax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`。LucaOneTasks README 对 SupKTax/GenusTax 的自然语言层级存在对调，因此解释结果时以实际 `label.txt` 为准。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="genustax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.metadata["predicted_label"])
```

- **示例结果（结构示意）：**

```text
(5,)
<概率最高的官方标签>
```

## 21) 物种分类

- **作用：** 预测 gene 序列对应的物种类别。
- **预测方法：** 调用 LucaOneTasks 官方 `SpeciesTax` 多分类 checkpoint 直接推理，返回 Top-k 物种标签和对应概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="speciestax", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`；可用 `top_k` 设置返回标签数。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("gene-1", "ATGGCC" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="speciestax",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个物种标签>,)
```

## 22) 原核蛋白亚细胞定位

- **作用：** 预测原核蛋白的亚细胞定位类别。
- **预测方法：** 调用 LucaOneTasks 官方 `ProtLoc` 多分类 checkpoint 直接推理，返回 Top-k 定位标签及概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_location", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="protein"` 的 `BiologicalSequence`。
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
        task="protein_location",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个亚细胞定位标签>,)
```

## 23) 蛋白稳定性

- **作用：** 根据蛋白序列输出官方 ProtStab 任务定义下的稳定性回归值。
- **预测方法：** 调用 LucaOneTasks 官方 `ProtStab` 回归 checkpoint 和 `predict_v1.py` 直接推理，不在 DNAKit 中拟合回归器。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_stability", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="protein"` 的 `BiologicalSequence`。
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

## 24) ncRNA 家族

- **作用：** 预测非编码 RNA 序列所属的家族。
- **预测方法：** 调用 LucaOneTasks 官方 `ncRNAFam` 多分类 checkpoint 直接推理，返回 Top-k 家族标签和概率。
- **API：** `dnakit.predictions.predict_sequence_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="ncrna_family", model_source_path=...)[必须], backend[可选])`
- **输入：** `sequence_type="gene"` 的 `BiologicalSequence`；DNA 或 RNA 字母均可，`U` 由官方任务流程处理。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [BiologicalSequence("ncrna-1", "AUGCUA" * 100, "gene")],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="ncrna_family",
        model_source_path="/opt/LucaOneTasks",
        top_k=5,
    ),
)
output = result.records[0].output
print(output.values.shape)
print(output.output_names)
```

- **示例结果（结构示意）：**

```text
(5,)
(<概率最高的 5 个 ncRNA 家族标签>,)
```

## 25) 甲型流感抗原关系

- **作用：** 预测两条甲型流感 gene 序列之间的抗原关系类别。
- **预测方法：** 调用 LucaOneTasks 官方 `InfA` 二分类 checkpoint 直接推理，返回正类概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="influenza_antigenicity", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；两条序列都必须是 `sequence_type="gene"`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "influenza-pair-1",
    BiologicalSequence("strain-a", "ATGGCC" * 100, "gene"),
    BiologicalSequence("strain-b", "ATGGCT" * 100, "gene"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="influenza_antigenicity",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```

## 26) 蛋白–蛋白互作

- **作用：** 预测两条蛋白序列是否具有蛋白–蛋白相互作用。
- **预测方法：** 调用 LucaOneTasks 官方 `PPI` 二分类 checkpoint 直接推理，返回互作概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="protein_interaction", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；两条序列都必须是 `sequence_type="protein"`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "protein-pair-1",
    BiologicalSequence("protein-a", "MAMAPRTEINSTRING" * 10, "protein"),
    BiologicalSequence("protein-b", "MKWVTFISLLFLFSSAYSR" * 8, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_interaction",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```

## 27) ncRNA–蛋白互作

- **作用：** 预测一条 ncRNA 与一条蛋白序列是否具有相互作用。
- **预测方法：** 调用 LucaOneTasks 官方 `ncRPI` 二分类 checkpoint 直接推理，返回互作概率和官方标签。
- **API：** `dnakit.predictions.predict_pair_properties(inputs[必须], config=PropertyPredictionConfig(model="lucaone", task="ncrna_protein_interaction", model_source_path=...)[必须], backend[可选])`
- **输入：** `BiologicalSequencePair`；第一条必须是 `sequence_type="gene"`，第二条必须是 `sequence_type="protein"`。
- **示例代码：**

```python
from dnakit.predictions import (
    BiologicalSequence,
    BiologicalSequencePair,
    PropertyPredictionConfig,
    predict_pair_properties,
)

pair = BiologicalSequencePair(
    "ncrna-protein-1",
    BiologicalSequence("ncrna-1", "AUGCUA" * 100, "gene"),
    BiologicalSequence("protein-1", "MKWVTFISLLFLFSSAYSR" * 8, "protein"),
)
result = predict_pair_properties(
    [pair],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="ncrna_protein_interaction",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.metadata["predicted_label"])
```

- **示例结果（标签由官方 checkpoint 决定）：**

```text
('probability',)
<官方预测标签>
```

## 28) H2A.Z 组蛋白变体区域判别

- **作用：** 判断输入序列是否对应 H2A.Z（H2AFZ）ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h2afz"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H2AFZ.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H2AFZ.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h2afz",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 29) H3K27ac 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K27ac ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k27ac"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K27ac.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K27ac.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k27ac",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 30) H3K27me3 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K27me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k27me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K27me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K27me3.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k27me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 31) H3K36me3 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K36me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k36me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K36me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K36me3.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k36me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 32) H3K4me1 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K4me1 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k4me1"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K4me1.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K4me1.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me1",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 33) H3K4me2 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K4me2 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k4me2"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K4me2.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K4me2.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me2",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 34) H3K4me3 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K4me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k4me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K4me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K4me3.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k4me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 35) H3K9ac 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K9ac ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k9ac"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K9ac.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K9ac.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k9ac",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 36) H3K9me3 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H3K9me3 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h3k9me3"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H3K9me3.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H3K9me3.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h3k9me3",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **本地示例 checkpoint 的 CUDA smoke 结果：**

```text
('0', '1')
[0.04569203 0.95430803]
1
```

## 37) H4K20me1 组蛋白修饰区域判别

- **作用：** 判断输入序列是否对应 H4K20me1 ChIP-seq 峰区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="h4k20me1"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`H4K20me1.ckpt`；默认放在 `./ckpt/enformer-benchmarks/H4K20me1.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 250)],
    task="h4k20me1",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 38) 增强子判别

- **作用：** 判断输入序列属于增强子还是非增强子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="enhancers"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`enhancers.ckpt`；默认放在 `./ckpt/enformer-benchmarks/enhancers.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="enhancers",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 39) 增强子类型分类

- **作用：** 在非增强子、组织特异增强子和组织不变增强子之间进行三分类。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 3 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="enhancers_types"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`enhancers_types.ckpt`；默认放在 `./ckpt/enformer-benchmarks/enhancers_types.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="enhancers_types",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>, <类别2概率>)
<概率最高的 checkpoint 标签>
```

## 40) 启动子判别

- **作用：** 判断输入序列属于启动子还是非启动子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="promoter_all"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`promoter_all.ckpt`；默认放在 `./ckpt/enformer-benchmarks/promoter_all.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_all",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 41) 非 TATA 启动子判别

- **作用：** 判断输入序列属于非 TATA 启动子还是负样本。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="promoter_no_tata"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`promoter_no_tata.ckpt`；默认放在 `./ckpt/enformer-benchmarks/promoter_no_tata.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_no_tata",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 42) TATA 启动子判别

- **作用：** 判断输入序列属于 TATA 启动子还是负样本。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="promoter_tata"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`promoter_tata.ckpt`；默认放在 `./ckpt/enformer-benchmarks/promoter_tata.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="promoter_tata",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 43) 剪接受体位点判别

- **作用：** 判断序列中央上下文是否包含剪接受体位点。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="splice_sites_acceptors"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`splice_sites_acceptors.ckpt`；默认放在 `./ckpt/enformer-benchmarks/splice_sites_acceptors.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_acceptors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 44) 剪接位点类型分类

- **作用：** 在无剪接位点、剪接受体和剪接供体之间进行三分类。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 3 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="splice_sites_all"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`splice_sites_all.ckpt`；默认放在 `./ckpt/enformer-benchmarks/splice_sites_all.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_all",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>, <类别2概率>)
<概率最高的 checkpoint 标签>
```

## 45) 剪接供体位点判别

- **作用：** 判断序列中央上下文是否包含剪接供体位点。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="splice_sites_donors"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`splice_sites_donors.ckpt`；默认放在 `./ckpt/enformer-benchmarks/splice_sites_donors.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 150)],
    task="splice_sites_donors",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 46) 编码区与基因间区判别

- **作用：** 判断输入来自蛋白编码区还是基因间区。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="demo_coding_vs_intergenomic_seqs"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`demo_coding_vs_intergenomic_seqs.ckpt`；默认放在 `./ckpt/enformer-benchmarks/demo_coding_vs_intergenomic_seqs.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 50)],
    task="demo_coding_vs_intergenomic_seqs",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 47) 人类与线虫序列判别

- **作用：** 判断输入基因组片段来自人类还是秀丽隐杆线虫。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="demo_human_or_worm"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`demo_human_or_worm.ckpt`；默认放在 `./ckpt/enformer-benchmarks/demo_human_or_worm.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 50)],
    task="demo_human_or_worm",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 48) 果蝇增强子判别

- **作用：** 判断输入是否属于 Stark 果蝇增强子集合。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="drosophila_enhancers_stark"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`drosophila_enhancers_stark.ckpt`；默认放在 `./ckpt/enformer-benchmarks/drosophila_enhancers_stark.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 500)],
    task="drosophila_enhancers_stark",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 49) 小鼠增强子判别

- **作用：** 判断输入是否属于 Ensembl 小鼠增强子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="dummy_mouse_enhancers_ensembl"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`dummy_mouse_enhancers_ensembl.ckpt`；默认放在 `./ckpt/enformer-benchmarks/dummy_mouse_enhancers_ensembl.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 500)],
    task="dummy_mouse_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 50) 人类 Cohn 增强子判别

- **作用：** 判断输入是否属于 Cohn 人类增强子集合。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_enhancers_cohn"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`human_enhancers_cohn.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_enhancers_cohn.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 125)],
    task="human_enhancers_cohn",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 51) 人类 Ensembl 增强子判别

- **作用：** 判断输入是否属于 Ensembl 人类增强子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_enhancers_ensembl"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`human_enhancers_ensembl.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_enhancers_ensembl.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 75)],
    task="human_enhancers_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```

## 52) 人类调控元件类型分类

- **作用：** 在人类增强子、启动子和开放染色质区域之间进行三分类。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 3 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_ensembl_regulatory"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`human_ensembl_regulatory.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_ensembl_regulatory.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 100)],
    task="human_ensembl_regulatory",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>, <类别2概率>)
<概率最高的 checkpoint 标签>
```

## 53) 人类非 TATA 启动子判别

- **作用：** 判断输入是否属于人类非 TATA 启动子。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_nontata_promoters"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`human_nontata_promoters.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_nontata_promoters.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", ("ACGT" * 63)[:251])],
    task="human_nontata_promoters",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **本地示例 checkpoint 的 CUDA smoke 结果：**

```text
('negative', 'positive')
[3.8891116e-08 1.0000000e+00]
positive
```

## 54) 人类开放染色质区域判别

- **作用：** 判断输入是否属于 Ensembl 人类开放染色质区域。
- **预测方法：** 加载该任务已经完整微调的 Enformer 主干与分类头，对最终层 embedding 做有效 bin 均值池化，再用 softmax 返回 2 个类别概率；运行时不训练或微调。
- **API：** `dnakit.predictions.predict_enformer_benchmark(inputs[必须], task="human_ocr_ensembl"[必须], checkpoint_dir或checkpoint_path[可选])`
- **输入：** `BiologicalSequence` 或 `DNARecord`；允许 A/C/G/T 和 IUPAC 模糊碱基，不会静默截断生物学序列。
- **Checkpoint：** 缺失时自动下载（见[权重缓存说明](#checkpoint-cache)）`human_ocr_ensembl.ckpt`；默认放在 `./ckpt/enformer-benchmarks/human_ocr_ensembl.ckpt`，也可使用自定义目录或精确路径。
- **示例代码：**

```python
from dnakit.predictions import BiologicalSequence, predict_enformer_benchmark

result = predict_enformer_benchmark(
    [BiologicalSequence("sequence-1", "ACGT" * 80)],
    task="human_ocr_ensembl",
    checkpoint_dir="/data/enformer-checkpoints",
    device="cuda",
)
output = result.records[0].output
print(output.output_names)
print(output.values)
print(output.metadata["predicted_label"])
```

- **示例结果（数值随序列与 checkpoint 而变化）：**

```text
(<类别0概率>, <类别1概率>)
<概率最高的 checkpoint 标签>
```
