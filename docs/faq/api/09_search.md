# 通用搜索：参数与补充说明

用法和示例见[功能教程](../../api/features/09_search.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 使用远程 BLAST 检索相似蛋白质序列 {#dnakit.databases.blast_protein_search-parameters}

`dnakit.databases.blast_protein_search(*, sequence, database=None, expect=None, hitlist_size=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `sequence` | string | 必填 | Protein sequence to search (minimum 10 amino acids) |
| `database` | string | None | Database to search. Options: &#x27;nr&#x27; (comprehensive but slow), &#x27;swissprot&#x27; (curated, faster), &#x27;pdb&#x27; (protein structures only) {&quot;default&quot;: &quot;nr&quot;} |
| `expect` | number | None | E-value threshold for reporting hits {&quot;default&quot;: 10.0} |
| `hitlist_size` | integer | None | Maximum number of hits to return {&quot;default&quot;: 50} |
| `verbose` | bool | True | 运行提示开关。 |
