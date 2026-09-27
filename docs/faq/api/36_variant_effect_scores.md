# 变异效应评分与实验数据：参数与补充说明

用法和示例见[功能教程](../../api/features/36_variant_effect_scores.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 查询整条蛋白质的 AlphaMissense 致病性分类 {#dnakit.databases.alpha_missense_get_protein_scores-parameters}

`dnakit.databases.alpha_missense_get_protein_scores(*, uniprot_id, max_residues=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 必填 | UniProt accession ID (e.g., &#x27;P00533&#x27; for EGFR, &#x27;P01116&#x27; for KRAS, &#x27;P04637&#x27; for TP53). |
| `max_residues` | integer &#124; null | None | Optional cap on the number of residues fetched. 0 / null / omitted = no cap (fetch every residue of the protein, recommended for whole-protein analysis). Positive integer = fetch only the first N residues (useful for fast diagnostic calls on very long proteins). {&quot;default&quot;: 0} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询蛋白质指定位点的替换分类与平均评分 {#dnakit.databases.alpha_missense_get_residue_scores-parameters}

`dnakit.databases.alpha_missense_get_residue_scores(*, uniprot_id, position, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 必填 | UniProt accession ID (e.g., &#x27;P00533&#x27; for EGFR) |
| `position` | integer | 必填 | Amino acid position in the protein (1-indexed) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询单个错义变异的 AlphaMissense 分类 {#dnakit.databases.alpha_missense_get_variant_score-parameters}

`dnakit.databases.alpha_missense_get_variant_score(*, uniprot_id, variant, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `uniprot_id` | string | 必填 | UniProt accession ID (e.g., &#x27;P00533&#x27; for EGFR) |
| `variant` | string | 必填 | Variant in protein notation: p.X123Y or X123Y where X is reference amino acid, 123 is position, Y is variant (e.g., &#x27;p.R123H&#x27;, &#x27;V600E&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因组指定位点的全部替换 CADD 分数 {#dnakit.databases.cadd_get_position_scores-parameters}

`dnakit.databases.cadd_get_position_scores(*, chrom, pos, version=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chrom` | string | 必填 | Chromosome (1-22, X, Y, MT) |
| `pos` | integer | 必填 | Genomic position (1-based) |
| `version` | string | None | CADD version and genome build {&quot;default&quot;: &quot;GRCh38-v1.7&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因组区间内的 CADD 分数 {#dnakit.databases.cadd_get_range_scores-parameters}

`dnakit.databases.cadd_get_range_scores(*, chrom, start, end, version=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chrom` | string | 必填 | Chromosome (1-22, X, Y, MT) |
| `start` | integer | 必填 | Start position (1-based, inclusive) |
| `end` | integer | 必填 | End position (1-based, inclusive). Max range: 100bp. |
| `version` | string | None | CADD version and genome build {&quot;default&quot;: &quot;GRCh38-v1.7&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 通过 MyVariant 查询已收录的 CADD 注释 {#dnakit.databases.cadd_get_scores-parameters}

`dnakit.databases.cadd_get_scores(*, variant, assembly='hg38', verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填 | 用于 MyVariant CADD 注释检索的变异查询条件。 |
| `assembly` | string | &#x27;hg38&#x27; | 参考基因组版本，例如 hg38。 {&quot;default&quot;: &quot;hg38&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询单个变异的 CADD 有害性分数 {#dnakit.databases.cadd_get_variant_score-parameters}

`dnakit.databases.cadd_get_variant_score(*, chrom, pos, ref, alt, version=None, include_annotations=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `chrom` | string | 必填 | Chromosome (1-22, X, Y, MT). Can include &#x27;chr&#x27; prefix. |
| `pos` | integer | 必填 | Genomic position (1-based) |
| `ref` | string | 必填 | Reference allele (e.g., &#x27;A&#x27;, &#x27;G&#x27;) |
| `alt` | string | 必填 | Alternate allele (e.g., &#x27;T&#x27;, &#x27;C&#x27;) |
| `version` | string | None | CADD version and genome build {&quot;default&quot;: &quot;GRCh38-v1.7&quot;, &quot;enum&quot;: &#91;&quot;GRCh38-v1.7&quot;, &quot;GRCh37-v1.7&quot;, &quot;GRCh38-v1.6&quot;, &quot;GRCh37-v1.6&quot;&#93;} |
| `include_annotations` | boolean | None | Include full annotation details in response {&quot;default&quot;: false} |
| `verbose` | bool | True | 运行提示开关。 |


## 查询基因的 EVE 评分覆盖情况 {#dnakit.databases.eve_get_gene_info-parameters}

`dnakit.databases.eve_get_gene_info(*, gene_symbol, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 必填 | Gene symbol (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;EGFR&#x27;) |
| `verbose` | bool | True | 运行提示开关。 |


## 查询变异的 EVE 致病性分数 {#dnakit.databases.eve_get_variant_score-parameters}

`dnakit.databases.eve_get_variant_score(*, variant=None, chrom=None, pos=None, ref=None, alt=None, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | None | Variant in HGVS format (e.g., &#x27;ENST00000269305.4:c.100G&gt;A&#x27;, &#x27;NM_000546.5:c.215C&gt;G&#x27;) |
| `chrom` | string | None | Chromosome (1-22, X, Y). Use with pos, ref, alt instead of variant. |
| `pos` | integer | None | Genomic position (GRCh38) |
| `ref` | string | None | Reference allele |
| `alt` | string | None | Alternate allele |
| `species` | string | None | Species (default: human) {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 将深度突变扫描数据整理为氨基酸效应矩阵 {#dnakit.databases.mave_db_get_effect_matrix-parameters}

`dnakit.databases.mave_db_get_effect_matrix(*, urn, uniprot_id=None, score_field=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `urn` | string | 必填 | MaveDB score set URN (e.g., &#x27;urn:mavedb:00000115-a-7&#x27;). Obtain from MaveDB_search_score_sets. |
| `uniprot_id` | string &#124; null | None | Optional UniProt accession (e.g., &#x27;P01116&#x27; for KRAS). If supplied, the tool fetches the canonical sequence and verifies the MaveDB position numbering against it via a landmark check; sets numbering_offset and surfaces a warning in numbering_check if the positions don&#x27;t match (silent off-by-N is the #1 silent bug in DMS analysis). |
| `score_field` | string &#124; null | None | Optional explicit column name to use as the variant score (e.g., &#x27;score&#x27;, &#x27;ddG&#x27;, &#x27;fitness&#x27;). If omitted, tries &#x27;score&#x27; → &#x27;ddG&#x27; → &#x27;fitness&#x27; → first numeric non-HGVS field in that order. The actual field used is reported in data.score_field_used. |
| `verbose` | bool | True | 运行提示开关。 |


## 汇总变异的多来源致病性预测分数 {#dnakit.databases.my_variant_get_pathogenicity_scores-parameters}

`dnakit.databases.my_variant_get_pathogenicity_scores(*, variant_id, fields=None, assembly=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填 | Variant ID: rsID (e.g., rs45478192) or HGVS genomic (e.g., chr16:g.23635348A&gt;C). Both formats work; an HGVS id must be expressed in the assembly given by the &#x27;assembly&#x27; parameter. |
| `fields` | string | None | Fields to return (pre-configured for pathogenicity scores) {&quot;default&quot;: &quot;dbnsfp.revel.score,cadd.phred,dbnsfp.alphamissense.score,dbnsfp.alphamissense.pred,dbnsfp.sift.score,dbnsfp.sift.pred,dbnsfp.polyphen2.hdiv.score,dbnsfp.polyphen2.hdiv.pred,dbnsfp.metarnn.score,dbnsfp.metarnn.pred,cadd.gerp.rs,dbnsfp.phylop.100way_vertebrate.rankscore,dbnsfp.phastcons.100way_vertebrate.rankscore,dbnsfp.vest4.score,dbnsfp.mutationtaster.pred,clinvar.rcv.clinical_significance,dbsnp.rsid&quot;} |
| `assembly` | string | None | Reference genome assembly for both the HGVS variant_id and the returned coordinates: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `verbose` | bool | True | 运行提示开关。 |
