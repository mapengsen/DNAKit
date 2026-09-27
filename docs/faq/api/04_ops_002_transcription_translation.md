# 转录与翻译：参数与补充说明

用法和示例见[功能教程](../../api/features/04_ops_002_transcription_translation.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 翻译三个正向阅读框并查找最长开放阅读框 {#dnakit.databases.dna_translate_reading_frames-parameters}

`dnakit.databases.dna_translate_reading_frames(*, operation=None, sequence, frame=None, genetic_code=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;translate_reading_frames&quot;&#93;, &quot;default&quot;: &quot;translate_reading_frames&quot;} |
| `sequence` | string | 必填 | DNA sequence (A, T, G, C only). Case insensitive. Spaces and newlines are stripped. |
| `frame` | string &#124; null | None | Reading frame to translate: &#x27;1&#x27;, &#x27;2&#x27;, &#x27;3&#x27;, or &#x27;all&#x27; (default). Frame 1 starts at position 0, frame 2 at position 1, frame 3 at position 2. {&quot;enum&quot;: &#91;&quot;1&quot;, &quot;2&quot;, &quot;3&quot;, &quot;all&quot;, null&#93;, &quot;default&quot;: &quot;all&quot;} |
| `genetic_code` | integer &#124; null | None | NCBI genetic code number. Only 1 (the standard code) is implemented; any other value is rejected with an error rather than translated with the standard table. Alternative code tables (2 = vertebrate mitochondrial, 3 = yeast mitochondrial, 11 = bacterial, etc.) are not supported by this tool. {&quot;enum&quot;: &#91;1, null&#93;, &quot;default&quot;: 1} |
| `verbose` | bool | True | 运行提示开关。 |
