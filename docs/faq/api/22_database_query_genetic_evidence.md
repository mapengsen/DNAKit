# 遗传病关联、剂量与可干预性：参数与补充说明

用法和示例见[功能教程](../../api/features/22_database_query_genetic_evidence.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 按基因查询 ClinGen 剂量敏感性审定结果 {#dnakit.databases.clin_gen_dosage_by_gene-parameters}

`dnakit.databases.clin_gen_dosage_by_gene(*, gene=None, gene_symbol=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Gene symbol to search for (e.g., &#x27;BRCA1&#x27;, &#x27;MECP2&#x27;, &#x27;TP53&#x27;). Case-insensitive. |
| `gene_symbol` | string | None | DNAKit optional alias for gene |
| `verbose` | bool | True | 运行提示开关。 |


## 按区间检索 ClinGen 剂量敏感性记录 {#dnakit.databases.clin_gen_dosage_region_search-parameters}

`dnakit.databases.clin_gen_dosage_region_search(*, chromosome, start, end, assembly=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chromosome` | string | 必填 | Chromosome (e.g., &#x27;17&#x27;, &#x27;X&#x27;). |
| `start` | integer | 必填 | Start position. |
| `end` | integer | 必填 | End position. |
| `assembly` | string | None | Genome assembly version. {&quot;enum&quot;: &#91;&quot;GRCh37&quot;, &quot;GRCh38&quot;&#93;, &quot;default&quot;: &quot;GRCh38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询成人遗传发现的临床可干预性证据 {#dnakit.databases.clin_gen_get_actionability_adult-parameters}

`dnakit.databases.clin_gen_get_actionability_adult(*, gene=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | 运行提示开关。 |


## 查询儿童遗传发现的临床可干预性证据 {#dnakit.databases.clin_gen_get_actionability_pediatric-parameters}

`dnakit.databases.clin_gen_get_actionability_pediatric(*, gene=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | 运行提示开关。 |


## 获取 ClinGen 剂量敏感性审定记录 {#dnakit.databases.clin_gen_get_dosage_sensitivity-parameters}

`dnakit.databases.clin_gen_get_dosage_sensitivity(*, gene=None, include_regions=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `include_regions` | boolean | None | Include recurrent CNV regions (default: false) {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |


## 获取 ClinGen 基因与疾病关联有效性审定记录 {#dnakit.databases.clin_gen_get_gene_validity-parameters}

`dnakit.databases.clin_gen_get_gene_validity(*, gene=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `verbose` | bool | True | 运行提示开关。 |


## 查询 ClinGen 变异致病性审定结果 {#dnakit.databases.clin_gen_get_variant_classifications-parameters}

`dnakit.databases.clin_gen_get_variant_classifications(*, gene=None, variant=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | None | Optional: Filter by gene symbol |
| `variant` | string | None | Optional: filter to one variant. Accepts a ClinGen canonical allele ID (&#x27;CA16020993&#x27;, with or without a &#x27;CAR:&#x27; prefix), a ClinVar VariationID (&#x27;370074&#x27;), full HGVS (&#x27;NM_000277.3:c.1315+1G&gt;T&#x27;), or a protein change (&#x27;p.Arg408Trp&#x27;). The form supplied selects the matching server-side filter. HGVS and protein-change strings match on substring, so a short fragment can return many rows -- combine with &#x27;gene&#x27; to scope them. |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因检索成人与儿童临床可干预性 {#dnakit.databases.clin_gen_search_actionability-parameters}

`dnakit.databases.clin_gen_search_actionability(*, gene, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填 | Gene symbol to search (e.g., BRCA1, MLH1, LDLR) |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因名称检索剂量敏感性评分 {#dnakit.databases.clin_gen_search_dosage_sensitivity-parameters}

`dnakit.databases.clin_gen_search_dosage_sensitivity(*, gene, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填 | Gene symbol to search (e.g., MECP2, PMP22, RAI1) |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因名称检索疾病关联有效性 {#dnakit.databases.clin_gen_search_gene_validity-parameters}

`dnakit.databases.clin_gen_search_gene_validity(*, gene, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene` | string | 必填 | Gene symbol to search (e.g., BRCA1, TP53, CFTR) |
| `verbose` | bool | True | 运行提示开关。 |
