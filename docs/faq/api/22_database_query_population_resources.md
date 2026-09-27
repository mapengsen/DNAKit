# 人群、样本与连锁查询：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_population_resources.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 列出 IGSR 与千人基因组数据集合 {#dnakit.databases.igsr_list_data_collections-parameters}

`dnakit.databases.igsr_list_data_collections(*, operation=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: list_data_collections). {&quot;enum&quot;: &#91;&quot;list_data_collections&quot;&#93;, &quot;default&quot;: &quot;list_data_collections&quot;} |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 50, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 检索千人基因组人群信息 {#dnakit.databases.igsr_search_populations-parameters}

`dnakit.databases.igsr_search_populations(*, operation=None, superpopulation=None, query=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_populations). {&quot;enum&quot;: &#91;&quot;search_populations&quot;&#93;, &quot;default&quot;: &quot;search_populations&quot;} |
| `superpopulation` | string &#124; null | None | Filter by superpopulation code: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian). Leave empty for all. |
| `query` | string &#124; null | None | Free-text search query to find populations by name or description (e.g., &#x27;Chinese&#x27;, &#x27;Yoruba&#x27;, &#x27;British&#x27;). |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 按人群或集合检索千人基因组样本 {#dnakit.databases.igsr_search_samples-parameters}

`dnakit.databases.igsr_search_samples(*, operation=None, population=None, data_collection=None, sample_name=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_samples). {&quot;enum&quot;: &#91;&quot;search_samples&quot;&#93;, &quot;default&quot;: &quot;search_samples&quot;} |
| `population` | string &#124; null | None | Filter by population code (e.g., &#x27;CHB&#x27;, &#x27;YRI&#x27;, &#x27;GBR&#x27;, &#x27;CEU&#x27;). Use IGSR_search_populations to find valid codes. |
| `data_collection` | string &#124; null | None | Filter by data collection title (e.g., &#x27;1000 Genomes 30x on GRCh38&#x27;, &#x27;1000 Genomes on GRCh38&#x27;). |
| `sample_name` | string &#124; null | None | Look up specific sample by name (e.g., &#x27;HG00096&#x27;, &#x27;NA12878&#x27;). |
| `limit` | integer | None | Maximum number of results to return. {&quot;default&quot;: 25, &quot;minimum&quot;: 1, &quot;maximum&quot;: 100} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 SNP 的群体特异性连锁不平衡代理变异 {#dnakit.databases.ldlink_get_proxies-parameters}

`dnakit.databases.ldlink_get_proxies(*, variant, population=None, r2_threshold=None, genome_build=None, limit=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填 | Query SNP rsID, e.g. &#x27;rs7903146&#x27;. |
| `population` | string &#124; null | None | 1000 Genomes population code (default &#x27;CEU&#x27;). Examples: &#x27;CEU&#x27; (European), &#x27;YRI&#x27; (African), &#x27;CHB&#x27; (East Asian), or super-populations &#x27;EUR&#x27;/&#x27;AFR&#x27;/&#x27;EAS&#x27;/&#x27;SAS&#x27;/&#x27;AMR&#x27;. |
| `r2_threshold` | number &#124; null | None | Minimum R2 to report (default 0.8; high LD). |
| `genome_build` | string &#124; null | None | &#x27;grch38&#x27; (default) or &#x27;grch37&#x27;. |
| `limit` | integer &#124; null | None | Max proxies to return (1-500, default 50). |
| `verbose` | bool | True | 运行提示开关。 |
