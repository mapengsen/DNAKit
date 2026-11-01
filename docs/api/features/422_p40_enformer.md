# Enformer 预测与指标

## 作用

为 P39/P40 整合批次提供明确的 DNA 操作接口，返回可序列化结果。

## 处理方法

函数调用其对应的上游操作；显式指定本地文件，新输出路径不得已存在。UCSC/黑名单下载仅在主动调用时联网。验证环境缺少 NUPACK 与 BLAST+，相关科学计算流程尚未实测。Enformer 使用与架构和适配器配置一致的本地权重。

## API

`dnakit.predictions.enformer_p40_functions`

| Function | Operation |
| --- | --- |
| `enformer_poisson_loss` | Compute upstream poisson_loss, including its epsilon and mean reduction. |
| `enformer_pearson_correlation` | Compute upstream per-sequence Pearson coefficient for batch x positions x tracks. |
| `enformer_channel_correlation` | Accumulate upstream MeanPearsonCorrCoefPerChannel over successive track batches. |
| `enformer_predict_tracks` | Predict the human or mouse tracks using a complete local Enformer checkpoint. |
| `enformer_predict_custom_tracks` | Predict user-defined tracks with a trained HeadAdapterWrapper checkpoint. |
| `enformer_predict_context_tracks` | Predict context-specific tracks using trained hypernetwork adapter weights. |
| `enformer_predict_context_attention` | Predict tracks with trained cross-attention context-token adapter weights. |

## 输入

DNA 字符串、JSON 字典/列表、本地基因组文件或模型权重。区间坐标、热力学参数见各函数 docstring；NUPACK 默认材料为 DNA。全部接口支持关键字参数 `verbose`。

## 输出

含完成状态、后端信息和操作结果的字典；文件操作返回输出路径，绘图返回图片路径及对应数据。缺少可选依赖时抛出 `BackendUnavailableError`。各批 manifest 单独记录未验收流程。

**示例代码：**

```python
from dnakit.predictions.enformer_p40_functions import enformer_poisson_loss
print(enformer_poisson_loss([1.0], [1.0], verbose=False)['values'])
```

**示例输出：**

```text
1.0
```
