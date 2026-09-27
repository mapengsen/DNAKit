# 序列编辑：参数与补充说明

用法和示例见[功能教程](../../api/features/04_ops_004_sequence_editing.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 模拟限制性内切酶消化与片段生成 {#dnakit.databases.dna_virtual_digest-parameters}

`dnakit.databases.dna_virtual_digest(*, operation=None, sequence, enzymes=None, circular=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation (optional; defaults to &#x27;virtual_digest&#x27; for this tool). {&quot;enum&quot;: &#91;&quot;virtual_digest&quot;&#93;, &quot;default&quot;: &quot;virtual_digest&quot;} |
| `sequence` | string | 必填 | DNA sequence (A, T, G, C, N only). Case insensitive. |
| `enzymes` | array &#124; null | None | List of enzyme names to use for digestion (e.g., &#91;&#x27;EcoRI&#x27;, &#x27;BamHI&#x27;&#93;). If null, uses all 24 NEB enzymes. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `circular` | boolean &#124; null | None | Treat sequence as circular DNA (e.g., plasmid). Default: false (linear). {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |
