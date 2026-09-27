# 功能目录与数量说明

[返回功能树](../api/features/function_tree.md)

> 更新日期：2026-09-13。按当前工作区 `src/dnakit` 的公开导出、模型/编码目录、CLI、Agent 目录和功能页核对，包含尚未提交的实现。按数据准备、处理、分析、评价、可视化和大语言模型工具归类，工程能力收录于 FAQ。编号沿用功能页或需求追踪编号，不表示连续的功能计数；这里按用户能力归类，不逐一列出配置类、结果类型、兼容别名或底层辅助函数。可选依赖、外部程序和服务要求见对应页面；列入树中表示已有接口与实现，不表示本机已安装或完成所有外部后端的实测。

## 核对口径与数量


本次对照了公开模块导出、`bioinformatics` 各子模块的公开函数/类方法、`tools/catalog.py`、CLI 命令及 Galaxy 包装器。功能树按任务归类，同一能力不按所属 Python 模块重复计数。

| 代码目录 | 本次核对结果 | 数量含义 |
| --- | --- | --- |
| `descriptors.DESCRIPTOR_NAMES_V1` | 240 | 固定 schema 的输出字段，包含不同统计量和理化表相关字段，不是 240 个独立 API |
| `bioinformatics.encoding.DNA_ENCODINGS` | 44 | 编码名称，各编码输出维数随参数/序列而异 |
| `representations.available_embedding_models()` | 11 | 表征模型，均提供全量微调与线性探测；需各模型运行环境与权重 |
| `predictions.available_machine_learning_models()` | 37 | 算法条目，分类/回归支持范围分别记录 |
| `predictions.available_prediction_tasks()` | 54 | 登记的直接预测任务，不含用户新训练的性质 |

生物信息学算法的公式、输入条件与限制以[补充说明](bioinformatics.md)和代码 docstring 为准。例如基因名称注释依赖参考比对；NB2 差异表达没有 DESeq2/edgeR 的完整收缩流程；约束优化采用有预算的搜索。这些能力均按实际实现范围列出。
