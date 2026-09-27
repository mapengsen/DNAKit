# 变异注释与标识转换

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

## 1) HGVS 变异 VEP 注释 {#1-hgvs-vep}

- **作用：** 通过 Ensembl VEP 返回 HGVS 变异的转录本后果和影响信息。
- **处理方法：** 将 HGVS 表达式提交 Ensembl VEP，保留转录本 consequence 记录。
- **API：** `dnakit.annotation.annotate_variant_vep(hgvs_notation[必须], species[可选])`
- **输入：** 一个 HGVS 表达式；默认物种为 `human`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import annotate_variant_vep

result = annotate_variant_vep("NM_000546.6:c.215C>G")
print(result.data)
```

- **示例结果：** Ensembl VEP 返回的结构化 consequence 记录。

## 2) rsID 的 VEP 注释 {#2-rsid-vep}

- **作用：** 使用 Ensembl VEP 注释一个 dbSNP rsID。
- **处理方法：** 将 rsID 提交 Ensembl VEP，返回对应后果注释。
- **API：** `dnakit.annotation.annotate_rsid_vep(rsid[必须], species[可选])`
- **输入：** 一个 rsID，例如 `rs28934578`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import annotate_rsid_vep

result = annotate_rsid_vep("rs28934578")
print(result.data)
```

- **示例结果：** 与该 rsID 对应的 VEP consequence 记录。

## 3) 变异标识转换 {#3}

- **作用：** 将 rsID、HGVS 等变异标识转换为 Ensembl 可提供的等价表示。
- **处理方法：** 查询 Ensembl 的标识重编码接口，获取等价表示。
- **API：** `dnakit.annotation.recode_variant(variant_id[必须], species[可选])`
- **输入：** 一个变异标识和可选物种。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import recode_variant

result = recode_variant("rs28934578")
print(result.data)
```

- **示例结果：** 对应的 HGVS、基因组位置和其他标识。

<!-- BEGIN GENERATED API: variant-annotation -->

<span id="variant-annotation"></span>

## 4) 查询 Ensembl 变异后果与完整注释 {#dnakit.databases.ensembl_var_get_variant_consequences}

- **作用：** 查询 Ensembl 变异后果与完整注释。
- **处理方法：** 通过 `EnsemblVariationExt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_var_get_variant_consequences(*, variant_id, species=None, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回字典，包含 `name`, `source`, `most_severe_consequence`, `ancestral_allele` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_var_get_variant_consequences(
    variant_id='rs429358',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "name": "<string>",
  "source": "<string>",
  "most_severe_consequence": "<string>",
  "ancestral_allele": "<string>",
  "minor_allele": "<string>",
  "MAF": "<number>",
  "synonyms": [
    "<string>"
  ],
  "mappings": [
    {
      "seq_region_name": "<string>",
      "start": "<integer>",
      "end": "<integer>",
      "allele_string": "<string>",
      "strand": "<integer>",
      "assembly_name": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_var_get_variant_consequences-parameters)

## 5) 使用 HGVS 表达式进行 VEP 变异注释 {#dnakit.databases.ensembl_vep_annotate_hgvs}

- **作用：** 使用 HGVS 表达式进行 VEP 变异注释。
- **处理方法：** 通过 `EnsemblVEP` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_vep_annotate_hgvs(*, hgvs_notation, species=None, verbose=True)`
- **输入：** 必须提供 `hgvs_notation`。
- **输出：** 返回字典，包含 `input`, `assembly_name`, `seq_region_name`, `start` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_vep_annotate_hgvs(
    hgvs_notation='BRAF:p.Val600Glu',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "input": "<string>",
  "assembly_name": "<string>",
  "seq_region_name": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "strand": "<integer>",
  "allele_string": "<string>",
  "most_severe_consequence": "<string>"
}
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_annotate_hgvs-parameters)

## 6) 使用 rsID 进行 VEP 变异注释 {#dnakit.databases.ensembl_vep_annotate_rsid}

- **作用：** 使用 rsID 进行 VEP 变异注释。
- **处理方法：** 通过 `EnsemblVEP` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_vep_annotate_rsid(*, variant_id=None, rsid=None, rs_id=None, species=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `input`, `assembly_name`, `seq_region_name`, `start` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_vep_annotate_rsid(
    variant_id='rs7903146',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "input": "<string>",
  "assembly_name": "<string>",
  "seq_region_name": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "allele_string": "<string>",
  "most_severe_consequence": "<string>",
  "transcript_consequences": [
    {
      "gene_symbol": "<string>",
      "variant_allele": "<string>",
      "consequence_terms": "<array>",
      "impact": "<string>",
      "amino_acids": "<string>",
      "sift_prediction": "<string>",
      "polyphen_prediction": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_annotate_rsid-parameters)

## 7) 按基因组位置进行 VEP 变异注释 {#dnakit.databases.ensembl_vep_region}

- **作用：** 按基因组位置进行 VEP 变异注释。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_vep_region(*, species, region, allele=None, AncestralAllele=None, Blosum62=None, CADD=None, Conservation=None, DisGeNET=None, GeneSplicer=None, GO=None, LoF=None, MaxEntScan=None, Phenotypes=None, SIFT=None, PolyPhen=None, verbose=True)`
- **输入：** 必须提供 `species`, `region`。
- **输出：** 返回记录列表；每条记录包含 `allele_string`, `consequence_terms`, `most_severe_consequence`, `transcript_consequences` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_vep_region(
    species='human',
    region='21:25891796-25891796',
    allele='T',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "allele_string": "<string>",
    "consequence_terms": [
      "<value>"
    ],
    "most_severe_consequence": "<string>",
    "transcript_consequences": [
      "<value>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_region-parameters)

## 8) 转换变异编号与 HGVS 表达形式 {#dnakit.databases.ensembl_vep_variant_recoder}

- **作用：** 转换变异编号与 HGVS 表达形式。
- **处理方法：** 通过 `EnsemblVEP` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_vep_variant_recoder(*, variant_id, species=None, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回记录列表；每条记录包含 `allele`, `input`, `id`, `hgvsg` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_vep_variant_recoder(
    variant_id='rs429358',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "allele": "<string>",
    "input": "<string>",
    "id": [
      "<string>"
    ],
    "hgvsg": [
      "<string>"
    ],
    "hgvsc": [
      "<string>"
    ],
    "hgvsp": [
      "<string>"
    ],
    "spdi": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.ensembl_vep_variant_recoder-parameters)

## 9) 获取 GRCh38 变异的综合功能注释 {#dnakit.databases.favor_annotate_variant}

- **作用：** 获取 GRCh38 变异的综合功能注释。
- **处理方法：** 通过 `FAVORVariantAnnotation` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.favor_annotate_variant(*, variant, verbose=True)`
- **输入：** 必须提供 `variant`。
- **输出：** 返回字典，包含 `found`, `variant`, `gene_consequence`, `allele_frequency` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.favor_annotate_variant(
    variant='19-44908822-C-T',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "found": "<boolean>",
  "variant": "<object>",
  "gene_consequence": "<object>",
  "allele_frequency": "<object>",
  "deleteriousness": "<object>",
  "conservation": "<object>",
  "clinical": "<object>",
  "regulatory": "<object>"
}
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.favor_annotate_variant-parameters)

## 10) 按 HGVS 或 rsID 查询综合变异注释 {#dnakit.databases.my_variant_get_variant_annotation}

- **作用：** 按 HGVS 或 rsID 查询综合变异注释。
- **处理方法：** 通过 `MyVariant` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.my_variant_get_variant_annotation(*, variant_id, fields=None, assembly=None, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回字典，包含 `_id`, `dbsnp`, `clinvar`, `cadd` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.my_variant_get_variant_annotation(
    variant_id='chr1:g.218631822G>A',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "_id": "<string>",
  "dbsnp": {
    "rsid": "<string>",
    "gene": "<object>",
    "alleles": [
      "<value>"
    ]
  },
  "clinvar": "<object>",
  "cadd": {
    "phred": "<number>",
    "consequence": "<string>"
  },
  "gnomad_genome": "<object>"
}
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.my_variant_get_variant_annotation-parameters)

## 11) 按基因、位置或临床意义检索变异 {#dnakit.databases.my_variant_query_variants}

- **作用：** 按基因、位置或临床意义检索变异。
- **处理方法：** 通过 `MyVariant` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.my_variant_query_variants(*, query=None, fields=None, size=None, assembly=None, variant_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `took`, `total`, `hits` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.my_variant_query_variants(
    query='rs58991260',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "took": "<integer>",
  "total": "<integer>",
  "hits": [
    {
      "_id": "<string>",
      "dbsnp": "<object>",
      "clinvar": "<object>",
      "cadd": "<object>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.my_variant_query_variants-parameters)

## 12) 校验并转换 HGVS 变异描述 {#dnakit.databases.variant_validator_validate_variant}

- **作用：** 校验并转换 HGVS 变异描述。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.variant_validator_validate_variant(*, genome_build, variant_description, select_transcripts, verbose=True)`
- **输入：** 必须提供 `genome_build`, `variant_description`, `select_transcripts`。
- **输出：** 返回字典，包含 `flag`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.variant_validator_validate_variant(
    genome_build='GRCh38',
    variant_description='NM_007294.4:c.5266dup',
    select_transcripts='NM_007294.4',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "flag": "<string>",
  "metadata": {
    "variantvalidator_version": "<string>",
    "variantvalidator_hgvs_version": "<string>",
    "vvdb_version": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/35_variant_annotation.md#dnakit.databases.variant_validator_validate_variant-parameters)

<!-- END GENERATED API: variant-annotation -->
