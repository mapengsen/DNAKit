# DNA 描述符：参数与补充说明

用法和示例见[功能教程](../../api/features/05_all_descriptors.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 统计指定碱基或氨基酸的数量与比例 {#dnakit.databases.sequence_count_residues-parameters}

`dnakit.databases.sequence_count_residues(*, operation=None, sequence=None, uniprot_id=None, residue, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;count_residues&quot;&#93;, &quot;default&quot;: &quot;count_residues&quot;} |
| `sequence` | string | None | Sequence string (DNA/RNA/protein). Optional if uniprot_id is provided. |
| `uniprot_id` | string | None | UniProt accession to fetch sequence from (e.g., &#x27;P00533&#x27; for EGFR). Optional if sequence is provided. |
| `residue` | string | 必填 | Single character residue to count (e.g., &#x27;C&#x27; for cysteine) |
| `verbose` | bool | True | 运行提示开关。 |


## 计算核酸序列的 GC 含量 {#dnakit.databases.sequence_gc_content-parameters}

`dnakit.databases.sequence_gc_content(*, operation=None, sequence, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;gc_content&quot;&#93;, &quot;default&quot;: &quot;gc_content&quot;} |
| `sequence` | string | 必填 | DNA or RNA sequence string |
| `verbose` | bool | True | 运行提示开关。 |


## 统计序列长度、组成与类型特征 {#dnakit.databases.sequence_stats-parameters}

`dnakit.databases.sequence_stats(*, operation=None, sequence=None, uniprot_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type {&quot;enum&quot;: &#91;&quot;stats&quot;&#93;, &quot;default&quot;: &quot;stats&quot;} |
| `sequence` | string | None | Sequence string. Optional if uniprot_id is provided. |
| `uniprot_id` | string | None | UniProt accession to fetch protein sequence from. Optional if sequence is provided. |
| `verbose` | bool | True | 运行提示开关。 |
