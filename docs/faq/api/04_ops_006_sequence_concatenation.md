# 序列拼接与组装：参数与补充说明

用法和示例见[功能教程](../../api/features/04_ops_006_sequence_concatenation.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 设计 Gibson 无缝组装重叠序列 {#dnakit.databases.dna_gibson_design-parameters}

`dnakit.databases.dna_gibson_design(*, operation=None, fragments, overlap_length=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;gibson_design&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;gibson_design&quot;&#93;, &quot;default&quot;: &quot;gibson_design&quot;} |
| `fragments` | array | 必填 | List of DNA fragment sequences (at least 2). Each must be longer than overlap_length. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `overlap_length` | integer &#124; null | None | Overlap length in bp for Gibson Assembly (default: 20 bp). Must be at least 1. {&quot;default&quot;: 20} |
| `verbose` | bool | True | 运行提示开关。 |


## 推演 Golden Gate 反应的组装产物 {#dnakit.databases.dna_golden_gate_assemble-parameters}

`dnakit.databases.dna_golden_gate_assemble(*, operation=None, fragments, enzyme=None, circular=None, labels=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;golden_gate_assemble&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;golden_gate_assemble&quot;&#93;, &quot;default&quot;: &quot;golden_gate_assemble&quot;} |
| `fragments` | array | 必填 | The DNA sequences combined in the reaction (plasmids or linear parts), A/T/G/C only. At least 2. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enzyme` | string &#124; null | None | Type IIS enzyme used in the reaction, e.g. &#x27;BsaI&#x27; (default), &#x27;BbsI&#x27;, &#x27;Esp3I&#x27;/&#x27;BsmBI&#x27;, &#x27;SapI&#x27;. {&quot;default&quot;: &quot;BsaI&quot;} |
| `circular` | boolean &#124; null | None | Whether the inputs are circular plasmids (default true). Set false for linear parts. {&quot;default&quot;: true} |
| `labels` | array &#124; null | None | Optional names for the inputs (e.g. plasmid names), echoed in the assembly order. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `verbose` | bool | True | 运行提示开关。 |


## 设计 Golden Gate 酶切位点与连接末端 {#dnakit.databases.dna_golden_gate_design-parameters}

`dnakit.databases.dna_golden_gate_design(*, operation=None, parts, enzyme=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;golden_gate_design&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;golden_gate_design&quot;&#93;, &quot;default&quot;: &quot;golden_gate_design&quot;} |
| `parts` | array | 必填 | List of DNA part sequences to assemble (at least 2). A, T, G, C only. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `enzyme` | string &#124; null | None | Type IIS restriction enzyme to use. Options: &#x27;BsaI&#x27; (default) or &#x27;BbsI&#x27;. {&quot;enum&quot;: &#91;&quot;BsaI&quot;, &quot;BbsI&quot;, null&#93;, &quot;default&quot;: &quot;BsaI&quot;} |
| `verbose` | bool | True | 运行提示开关。 |
