# 变异注释与标识转换：参数与补充说明

用法和示例见[功能教程](../../api/features/35_variant_annotation.md)。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

各功能按作用、API、输入、输出、示例代码和结果说明介绍。运行前按说明准备凭据、本地文件与可选依赖；远程结果随数据库版本与查询时间变化。

## 查询 Ensembl 变异后果与完整注释 {#dnakit.databases.ensembl_var_get_variant_consequences-parameters}

`dnakit.databases.ensembl_var_get_variant_consequences(*, variant_id, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填 | dbSNP rsID of the variant. Examples: &#x27;rs429358&#x27; (APOE), &#x27;rs7903146&#x27; (TCF7L2), &#x27;rs1042779&#x27;. |
| `species` | string | None | Species name. Default &#x27;human&#x27;. Examples: &#x27;human&#x27;, &#x27;homo_sapiens&#x27;. |
| `verbose` | bool | True | 运行提示开关。 |


## 使用 HGVS 表达式进行 VEP 变异注释 {#dnakit.databases.ensembl_vep_annotate_hgvs-parameters}

`dnakit.databases.ensembl_vep_annotate_hgvs(*, hgvs_notation, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `hgvs_notation` | string | 必填 | Variant in HGVS notation. Supports protein (p.), coding DNA (c.), and genomic (g.) notation. Examples: &#x27;BRAF:p.Val600Glu&#x27;, &#x27;ENST00000366667:c.803C&gt;T&#x27;, &#x27;NC_000007.14:g.140753336A&gt;T&#x27;. |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. Other options: &#x27;mouse&#x27;, &#x27;rat&#x27;, &#x27;zebrafish&#x27;, etc. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 使用 rsID 进行 VEP 变异注释 {#dnakit.databases.ensembl_vep_annotate_rsid-parameters}

`dnakit.databases.ensembl_vep_annotate_rsid(*, variant_id=None, rsid=None, rs_id=None, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | None | dbSNP rs identifier. Examples: &#x27;rs7903146&#x27; (TCF7L2, type 2 diabetes), &#x27;rs429358&#x27; (APOE, Alzheimer&#x27;s), &#x27;rs1801133&#x27; (MTHFR C677T). Aliases: rsid, rs_id. |
| `rsid` | string | None | Alias for variant_id. dbSNP rs identifier (e.g., &#x27;rs7903146&#x27;). |
| `rs_id` | string | None | Alias for variant_id. dbSNP rs identifier (e.g., &#x27;rs7903146&#x27;). |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因组位置进行 VEP 变异注释 {#dnakit.databases.ensembl_vep_region-parameters}

`dnakit.databases.ensembl_vep_region(*, species, region, allele=None, AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, GeneSplicer=None, GO=None, LoF=None, MaxEntScan=None, Phenotypes=None, SIFT=None, PolyPhen=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `species` | string | 必填 | Species (e.g., &#x27;human&#x27;, &#x27;homo_sapiens&#x27;) |
| `region` | string | 必填 | Genomic region with variant (e.g., &#x27;21:25891796-25891796&#x27;) |
| `allele` | string | None | VEP 区域变异的等位基因参数，遵循 Ensembl VEP 格式。 |
| `AncestralAllele` | integer | None | Include ancestral allele information {&quot;default&quot;: 0} |
| `Blosum62` | integer | None | Include BLOSUM62 protein conservation score {&quot;default&quot;: 0} |
| `CADD` | integer | None | Include CADD scores {&quot;default&quot;: 0} |
| `Conservation` | integer | None | Include conservation scores {&quot;default&quot;: 0} |
| `DisGeNET` | integer | None | Include DisGeNET disease associations {&quot;default&quot;: 0} |
| `GeneSplicer` | integer | None | Include GeneSplicer predictions {&quot;default&quot;: 0} |
| `GO` | integer | None | Include Gene Ontology terms {&quot;default&quot;: 0} |
| `LoF` | integer | None | Include Loss-of-Function predictions {&quot;default&quot;: 0} |
| `MaxEntScan` | integer | None | Include MaxEntScan splice predictions {&quot;default&quot;: 0} |
| `Phenotypes` | integer | None | Include phenotype annotations {&quot;default&quot;: 0} |
| `SIFT` | string | None | SIFT prediction: p=prediction, s=score, b=both {&quot;enum&quot;: &#91;&quot;p&quot;, &quot;s&quot;, &quot;b&quot;&#93;} |
| `PolyPhen` | string | None | PolyPhen prediction: p=prediction, s=score, b=both {&quot;enum&quot;: &#91;&quot;p&quot;, &quot;s&quot;, &quot;b&quot;&#93;} |
| `verbose` | bool | True | 运行提示开关。 |


## 转换变异编号与 HGVS 表达形式 {#dnakit.databases.ensembl_vep_variant_recoder-parameters}

`dnakit.databases.ensembl_vep_variant_recoder(*, variant_id, species=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填 | Variant identifier to recode. Accepts: rsID (e.g., &#x27;rs429358&#x27;), HGVS notation (e.g., &#x27;AGT:c.803T&gt;C&#x27;), SPDI (e.g., &#x27;NC_000019.10:44908683:T:C&#x27;). |
| `species` | string | None | Species name. Default: &#x27;human&#x27;. {&quot;default&quot;: &quot;human&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 获取 GRCh38 变异的综合功能注释 {#dnakit.databases.favor_annotate_variant-parameters}

`dnakit.databases.favor_annotate_variant(*, variant, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant` | string | 必填 | GRCh38 SNV/indel as chr-pos-ref-alt. Accepts &#x27;19-44908822-C-T&#x27;, &#x27;chr19:44908822:C:T&#x27;, etc. |
| `verbose` | bool | True | 运行提示开关。 |


## 按 HGVS 或 rsID 查询综合变异注释 {#dnakit.databases.my_variant_get_variant_annotation-parameters}

`dnakit.databases.my_variant_get_variant_annotation(*, variant_id, fields=None, assembly=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `variant_id` | string | 必填 | HGVS variant ID (e.g., &#x27;chr1:g.2186318G&gt;A&#x27;) or rsID (e.g., &#x27;rs123&#x27;). An HGVS id must be expressed in the assembly given by the &#x27;assembly&#x27; parameter. |
| `fields` | string | None | Fields to return. Useful fields: dbsnp, clinvar, cadd, gnomad_genome, gnomad_exome, dbnsfp (functional predictions). {&quot;default&quot;: &quot;dbsnp,clinvar,cadd,gnomad_genome,dbnsfp&quot;} |
| `assembly` | string | None | Reference genome assembly for both the HGVS variant_id and the returned coordinates: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `verbose` | bool | True | 运行提示开关。 |


## 按基因、位置或临床意义检索变异 {#dnakit.databases.my_variant_query_variants-parameters}

`dnakit.databases.my_variant_query_variants(*, query=None, fields=None, size=None, assembly=None, variant_id=None, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Examples: &#x27;rs58991260&#x27; (rsID), &#x27;chr7:g.55249071G&gt;A&#x27; (HGVS), &#x27;clinvar.gene.symbol:BRCA1&#x27; (gene in ClinVar), &#x27;cadd.phred:&gt;30&#x27; (high CADD score). |
| `fields` | string | None | Comma-separated fields to return. Common: dbsnp, clinvar, cadd, gnomad_genome, dbnsfp. {&quot;default&quot;: &quot;dbsnp.rsid,clinvar.rcv.clinical_significance,cadd.phred,gnomad_genome.af.af&quot;} |
| `size` | integer | None | Maximum number of results (1-100). {&quot;default&quot;: 10} |
| `assembly` | string | None | Reference genome assembly for the coordinates in the query and in the results: &#x27;hg19&#x27; (GRCh37, the default) or &#x27;hg38&#x27; (GRCh38). {&quot;enum&quot;: &#91;&quot;hg19&quot;, &quot;hg38&quot;&#93;, &quot;default&quot;: &quot;hg19&quot;} |
| `variant_id` | string | None | DNAKit optional alias for query |
| `verbose` | bool | True | 运行提示开关。 |


## 校验并转换 HGVS 变异描述 {#dnakit.databases.variant_validator_validate_variant-parameters}

`dnakit.databases.variant_validator_validate_variant(*, genome_build, variant_description, select_transcripts, verbose=True)`

| 参数 | 类型 | 必填 / 默认值 | 说明与约束 |
| --- | --- | --- | --- |
| `genome_build` | string | 必填 | Reference genome assembly: &#x27;GRCh37&#x27; (hg19) or &#x27;GRCh38&#x27; (hg38) |
| `variant_description` | string | 必填 | HGVS variant description (e.g., &#x27;NM_007294.4:c.5266dup&#x27; for BRCA1 c.5266dupC, &#x27;NM_000179.3:c.943C&gt;T&#x27;, &#x27;chr17:g.43092919del&#x27;) |
| `select_transcripts` | string | 必填 | Transcript to validate against (e.g., &#x27;NM_007294.4&#x27;). Use &#x27;all&#x27; to get all transcripts for a genomic variant. |
| `verbose` | bool | True | 运行提示开关。 |
