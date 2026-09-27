# 序列补全：参数与补充说明

用法和示例见[功能教程](../../api/features/04_ops_001_sequence_direction.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 计算 DNA 序列的反向互补 {#dnakit.databases.sequence_reverse_complement-parameters}

`dnakit.databases.sequence_reverse_complement(*, operation=None, sequence, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;reverse_complement&quot;&#93;, &quot;default&quot;: &quot;reverse_complement&quot;} |
| `sequence` | string | 必填 | DNA sequence (A/T/C/G/N only) |
| `verbose` | bool | True | 运行提示开关。 |
