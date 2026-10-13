# G7 数据库函数：gget、MyGene、MyVariant 与 biothings_client

本页记录 `databases` 区域第 G7 组的包装函数。每个函数都从 `dnakit.databases` 导出，例如 `from dnakit.databases import gget_cosmic_local`。

- 本地函数：`gget_cosmic_local` 只读取本地 COSMIC 文件，不需要网络或账号。
- 网络函数：其余函数访问公共数据库。默认单元测试使用替身客户端，不发起真实请求；下方示例的输出为真实运行结果，因此需要能访问对应服务。
- 示例中的输入文件以 `示例文件` 标签给出；其余示例直接运行。

## 已有包装与未包装的函数

下列函数已经由 `gene_backends.py` 包装，因此在本组中标记为重复，不再单独包装：gget 的 `search`、`info`、`seq`、`mutate`；MyGene 的 `query`、`querymany`、`getgene`；MyVariant 的 `getvariants`、`query`、`format_hgvs`、`get_hgvs_from_vcf`；biothings_client 的 HGVS 辅助函数。`findgenes` 是已弃用别名。

其他未包装的函数分为：需要外部程序（`diamond`、`muscle`、`elm` 的 DIAMOND 路径）；需要大规模下载（`ref`、`virus`）；需要 GPU（`alphafold`）；需要额外包（`cellxgene` 需要 `cellxgene_census`）；绘图或 CLI、缓存与配置辅助函数；以及蛋白质层面的 `g2p`、`pdb`。完整状态见 G7 台账。

## 1) `gget_cosmic_local`

- **作用：** 在本地 COSMIC 突变表中检索基因名、突变编号、CDS 或 AA 变化等记录，不下载数据、不需要登录。
- **处理方法：** 调用 `gget.cosmic` 的本地分支。`cosmic_project` 缺省时与 gget 一致：文件名含 `CancerMutationCensus_AllData` 视为 `cancer`，否则为 `other`。无匹配返回空列表；列名中的空格替换为下划线。
- **API：** `gget_cosmic_local(searchterm[必须], cosmic_tsv_path[必须], cosmic_project[可选], limit[可选], verbose[可选])`
- **输入：** 检索词字符串；本地 COSMIC TSV 路径；`cosmic_project` 取 cancer、cancer_example、census、resistance、cell_line、genome_screen、targeted_screen 或 other；`limit` 为 1 至 10000。
- **输出：** 字典，含 `count`、`records`（按文件顺序排列）、`cosmic_project` 与 `provenance`。COSMIC 商业使用需遵守其许可条款。

**示例文件：** `CancerMutationCensus_AllData_mini.tsv`

```text
GENE_NAME	ACCESSION_NUMBER	LEGACY_MUTATION_ID	Mutation CDS	Mutation AA	GENOMIC_MUTATION_ID	Primary site
EGFR	NM_005228	COSM6223	c.2573T>G	p.L858R	COSV55365	lung
BRAF	NM_004333	COSM476	c.1799T>A	p.V600E	COSV55410	skin
EGFR	NM_005228	COSM6224	c.2236_2250del	p.E746_A750del	COSV55366	lung
```

**示例代码：**

```python
from dnakit.databases import gget_cosmic_local

result = gget_cosmic_local("EGFR", "CancerMutationCensus_AllData_mini.tsv", verbose=False)
print(result["count"], [row["Mutation_AA"] for row in result["records"]])
```

**示例输出：**

```text
2 ['p.L858R', 'p.E746_A750del']
```

## 2) `gget_blast`

- **作用：** 把一条核酸或蛋白序列提交到 NCBI BLAST（经 gget），返回命中列表。
- **处理方法：** 提交前统一大写并校验字母表、`program`、`database`、`limit` 和 `expect`。远程作业可能需要数分钟，因此示例只演示参数校验；网络路径在单元测试中用替身替代。
- **API：** `gget_blast(sequence[必须], program[可选], database[可选], limit[可选], expect[可选], verbose[可选])`
- **输入：** 序列字符串，最长 20000 个字符；`program` 取 default、blastn、blastp、blastx、tblastn 或 tblastx；`expect` 在 1e-30 至 1000 之间。
- **输出：** 字典，含 `query_length`、`program`、`database`、`expect`、`count` 与 `hits`。

**示例代码：**

```python
from dnakit.databases import gget_blast
from dnakit.exceptions import ConfigurationError

try:
    gget_blast("ACGTACGT", program="blastz", verbose=False)
except ConfigurationError as error:
    print(error)
```

**示例输出：**

```text
[CONFIGURATION_ERROR] program must be one of default, blastn, blastp, blastx, tblastn, tblastx; got 'blastz'.
```

## 3) `gget_blat`

- **作用：** 用 UCSC BLAT 在指定组装中定位 DNA 或蛋白序列。
- **处理方法：** `seqtype` 缺省时由 gget 按字母判断核酸或蛋白。包装器只校验参数并把结果转为普通字典；示例序列取自 hg38 chr17:7668402-7668461，因此第一条命中应与之完全一致。
- **API：** `gget_blat(sequence[必须], seqtype[可选], assembly[可选], verbose[可选])`
- **输入：** 序列字符串；`seqtype` 取 default、DNA、protein、translated RNA 或 translated DNA；`assembly` 为 UCSC 组装名，如 human。
- **输出：** 字典，含 `query_length`、`seqtype`、`assembly`、`count` 与 `hits`。每条命中包含 `chromosome`、`start`、`end`、`strand` 和 `%_matched` 等字段。

**示例代码：**

```python
from dnakit.databases import gget_blat

sequence = "CACCCCTCAGACACACAGGTGGCAGCAAAGTTTTATTGTAAAATAAGAGATCGATATAAA"
result = gget_blat(sequence, seqtype="DNA", verbose=False)
top = result["hits"][0]
print(result["query_length"], result["count"], top["chromosome"], top["start"], top["end"], top["%_matched"])
```

**示例输出：**

```text
60 4 chr17 7668402 7668461 100.0
```

## 4) `gget_opentargets`

- **作用：** 按 Ensembl 基因 ID 查询 OpenTargets 中与该基因相关的疾病、药物、表达、互作或药物基因组学记录。
- **处理方法：** 只接受 `ENS` 开头、以 `G` 表示基因的 Ensembl ID，`resource` 必须在支持列表中；调用 gget 的 OpenTargets GraphQL 接口。
- **API：** `gget_opentargets(ensembl_id[必须], resource[可选], limit[可选], verbose[可选])`
- **输入：** Ensembl 基因 ID，例如 ENSG00000157764（BRAF）；`resource` 取 depmap、diseases、drugs、expression、interactions、pharmacogenetics 或 tractability。
- **输出：** 字典，含 `ensembl_id`、`resource`、`count` 与 `records`。记录字段由 OpenTargets 返回，并以点号展开。

**示例代码：**

```python
from dnakit.databases import gget_opentargets

result = gget_opentargets("ENSG00000157764", resource="drugs", limit=2, verbose=False)
print(result["count"], [row["drug.name"] for row in result["records"]])
```

**示例输出：**

```text
2 ['BELVARAFENIB', 'PLIXORAFENIB']
```

## 5) `gget_archs4`

- **作用：** 在 ARCHS4 中查找与某基因共表达最相关的基因，或查看该基因的组织表达图谱。
- **处理方法：** `which` 取 correlation（共表达）或 tissue（组织表达）；`species` 取 human 或 mouse。`ensembl=True` 时要求输入 Ensembl 基因 ID。
- **API：** `gget_archs4(gene[必须], ensembl[可选], which[可选], gene_count[可选], species[可选], verbose[可选])`
- **输入：** 基因符号或 Ensembl ID；`gene_count` 为 1 至 1000。
- **输出：** 字典，含 `gene`、`which`、`species`、`count` 与 `records`。

**示例代码：**

```python
from dnakit.databases import gget_archs4

result = gget_archs4("SOX2", which="correlation", gene_count=3, verbose=False)
print(result["count"], [row["gene_symbol"] for row in result["records"]])
```

**示例输出：**

```text
3 ['SOX21', 'PTPRZ1', 'SALL3']
```

## 6) `gget_bgee`

- **作用：** 从 Bgee 获取一个基因的直系同源关系，或多个基因的表达调用结果。
- **处理方法：** `record_type='orthologs'` 只接受一个 Ensembl 基因 ID；`record_type='expression'` 接受最多 100 个 ID。
- **API：** `gget_bgee(gene_id[必须], record_type[可选], verbose[可选])`
- **输入：** 单个 Ensembl 基因 ID，或 ID 列表；`record_type` 取 orthologs 或 expression。
- **输出：** 字典，含 `gene_ids`、`record_type`、`count` 与 `records`。

**示例代码：**

```python
from dnakit.databases import gget_bgee

result = gget_bgee("ENSG00000139618", record_type="orthologs", verbose=False)
print(result["count"], result["records"][0]["gene_name"])
```

**示例输出：**

```text
131 BRCA2
```

## 7) `gget_enrichr`

- **作用：** 用 Enrichr 对基因列表做基因集富集分析，返回各通路或功能条目的统计结果。
- **处理方法：** 基因列表最多 500 个，去重与顺序保持不变；可选背景基因列表最多 50000 个。服务端排序的第一条即 `records[0]`。
- **API：** `gget_enrichr(genes[必须], database[必须], species[可选], background_list[可选], verbose[可选])`
- **输入：** 基因符号列表；`database` 为 Enrichr 库名，如 `KEGG_2021_Human`；`species` 取 human 或 mouse。
- **输出：** 字典，含 `gene_count`、`database`、`species`、`count` 与 `records`。

**示例代码：**

```python
from dnakit.databases import gget_enrichr

result = gget_enrichr(["TP53", "MDM2", "CDKN1A", "BAX"], database="KEGG_2021_Human", verbose=False)
print(result["gene_count"], result["count"], result["records"][0]["path_name"])
```

**示例输出：**

```text
4 71 Melanoma
```

## 8) `gget_specificity`

- **作用：** 从 8cubeDB 获取基因层面的组织或细胞特异性统计。
- **处理方法：** 基因列表校验后直接传给 gget（8cubeDB 接口）；返回记录保持服务端顺序。
- **API：** `gget_specificity(gene_list[必须], verbose[可选])`
- **输入：** 基因符号或 ID 列表，最多 500 个。
- **输出：** 字典，含 `gene_count`、`count` 与 `records`。

**示例代码：**

```python
from dnakit.databases import gget_specificity

result = gget_specificity(["SOX2", "PAX6"], verbose=False)
print(result["gene_count"], result["count"], result["records"][0]["gene_name"])
```

**示例输出：**

```text
2 124 Pax6
```

## 9) `gget_gene_expression`

- **作用：** 从 8cubeDB 获取归一化表达值，并按分析层级与分组设计筛选。
- **处理方法：** `analysis_level` 与 `analysis_type` 原样传给 8cubeDB（例如 Across_tissues 与 Sex）；两者只做字符校验。
- **API：** `gget_gene_expression(gene_list[必须], analysis_level[必须], analysis_type[必须], verbose[可选])`
- **输入：** 基因列表，最多 500 个；分析层级和分组设计字符串。
- **输出：** 字典，含 `gene_count`、`analysis_level`、`analysis_type`、`count` 与 `records`。

**示例代码：**

```python
from dnakit.databases import gget_gene_expression

result = gget_gene_expression(
    ["SOX2", "PAX6"], analysis_level="Across_tissues", analysis_type="Sex", verbose=False
)
print(result["count"], [row["gene_name"] for row in result["records"]])
```

**示例输出：**

```text
2 ['Pax6', 'Sox2']
```

## 10) `gget_psi_block`

- **作用：** 从 8cubeDB 获取 psi-block 特异性分数，用于比较基因在各分组中的特异性。
- **处理方法：** 与 `gget_gene_expression` 相同的校验规则，调用 8cubeDB 的 psi-block 接口。
- **API：** `gget_psi_block(gene_list[必须], analysis_level[必须], analysis_type[必须], verbose[可选])`
- **输入：** 基因列表，最多 500 个；分析层级和分组设计字符串。
- **输出：** 字典，含 `gene_count`、`analysis_level`、`analysis_type`、`count` 与 `records`。

**示例代码：**

```python
from dnakit.databases import gget_psi_block

result = gget_psi_block(
    ["SOX2", "PAX6"], analysis_level="Across_tissues", analysis_type="Sex", verbose=False
)
print(result["count"], [row["gene_name"] for row in result["records"]])
```

**示例输出：**

```text
2 ['Pax6', 'Sox2']
```

## 11) `mygene_getgenes`

- **作用：** 用 Entrez 或 Ensembl ID 在 MyGene.info 中批量获取基因记录。
- **处理方法：** 调用 `mygene.MyGeneInfo().getgenes`。结果按输入顺序排列；MyGene 无法识别的 ID 单独列在 `not_found` 中，不静默丢弃。
- **API：** `mygene_getgenes(gene_ids[必须], fields[可选], species[可选], verbose[可选])`
- **输入：** Entrez 或 Ensembl ID 列表，最多 1000 个；`fields` 为逗号分隔字段名；`species` 为 human 或 NCBI taxon ID。
- **输出：** 字典，含 `requested`、`found`、`not_found`、`records` 与 `provenance`。

**示例代码：**

```python
from dnakit.databases import mygene_getgenes

result = mygene_getgenes(["1017", "999999999", "7157"], fields="symbol,name", verbose=False)
print(result["found"], [row["symbol"] for row in result["records"]], result["not_found"])
```

**示例输出：**

```text
2 ['CDK2', 'TP53'] ['999999999']
```

## 12) `myvariant_querymany`

- **作用：** 用 HGVS 变异 ID 或 rsID 在 MyVariant.info 中批量查询变异注释。
- **处理方法：** 调用 `myvariant.MyVariantInfo().querymany`。同一输入可能返回多条命中，输出按输入顺序排列；`matched_terms` 统计至少有一条命中的输入项，`not_found` 列出没有命中的输入项。
- **API：** `myvariant_querymany(query_terms[必须], scopes[可选], fields[可选], assembly[可选], verbose[可选])`
- **输入：** HGVS ID 或 rsID 列表，最多 1000 个；`assembly` 为 hg19 或 hg38。
- **输出：** 字典，含 `requested`、`matched_terms`、`hit_count`、`not_found`、`records` 与 `provenance`。

**示例代码：**

```python
from dnakit.databases import myvariant_querymany

result = myvariant_querymany(["rs58991260"], scopes="dbsnp.rsid", fields="dbsnp.rsid", assembly="hg38", verbose=False)
print(result["matched_terms"], result["records"][0]["_id"], result["not_found"])
```

**示例输出：**

```text
1 chr1:g.218458480G>A []
```

## 13) `biothings_client_get_homologs`

- **作用：** 通过 biothings_client 的 `utils.gene.get_homologs` 查找基因在其他物种中的同源基因。
- **处理方法：** 先读取 homologene 注释，再批量取回同源基因记录。真实的连接逻辑由 biothings_client 完成；物种过滤在其结果上生效。
- **API：** `biothings_client_get_homologs(gene_list[必须], fields[可选], species[可选], verbose[可选])`
- **输入：** Entrez 或 Ensembl 基因 ID 列表，最多 500 个；`species` 为 `all` 或逗号分隔的 NCBI taxon ID。
- **输出：** 字典，含 `requested`、`results` 与 `provenance`。每个基因列出按 taxon ID 分组的同源基因记录。

**示例代码：**

```python
from dnakit.databases import biothings_client_get_homologs

result = biothings_client_get_homologs(["7157"], fields="symbol", species="10090", verbose=False)
homologs = result["results"][0]["homologs"]
print(result["requested"], [(h["taxid"], g["symbol"]) for h in homologs for g in h["genes"]])
```

**示例输出：**

```text
1 [(10090, 'Trp53')]
```
