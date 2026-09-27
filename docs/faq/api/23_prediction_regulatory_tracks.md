# 表达、调控与接触图预测：参数与补充说明

用法和示例见[功能教程](../../api/features/23_prediction_regulatory_tracks.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

这些入口调用已有模型或预测服务。模型权重、SDK、服务凭据和使用条件继续由原提供方决定；DNAKit 没有重新训练这些模型。

## 预测基因组区间的多模态调控轨道 {#dnakit.databases.alpha_genome_predict_interval-parameters}

`dnakit.databases.alpha_genome_predict_interval(*, chromosome, start, end, output_types=None, ontology_terms=None, organism=None, sequence_length=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chromosome` | string | 必填 | Chromosome, e.g. &#x27;chr19&#x27;. |
| `start` | integer | 必填 | Interval start (0-based). |
| `end` | integer | 必填 | Interval end. |
| `output_types` | array &#124; null | None | Modalities to predict, e.g. &#91;&#x27;RNA_SEQ&#x27;,&#x27;ATAC&#x27;&#93; (default &#91;&#x27;RNA_SEQ&#x27;&#93;). {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `ontology_terms` | array &#124; null | None | Optional tissue/cell ontology terms (e.g. &#91;&#x27;UBERON:0001114&#x27;&#93; = liver) to restrict tracks. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `organism` | string &#124; null | None | &#x27;human&#x27; (default) or &#x27;mouse&#x27;. |
| `sequence_length` | string &#124; null | None | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `verbose` | bool | True | 运行提示开关。 |

运行前配置：`ALPHA_GENOME_API_KEY`。
