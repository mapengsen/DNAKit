# 结构变异查询

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: structural-variants -->

<span id="structural-variants"></span>

## 1) 查询 Ensembl 区间内的结构变异 {#dnakit.databases.ensembl_get_structural_variants}

- **作用：** 查询 Ensembl 区间内的结构变异。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_get_structural_variants(*, species, region=None, feature=None, chrom=None, start=None, end=None, verbose=True)`
- **输入：** 必须提供 `species`。
- **输出：** 返回记录列表；每条记录包含 `id`, `seq_region_name`, `start`, `end` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_get_structural_variants(
    species='human',
    region='22:27000000-27050000',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "id": "<string>",
    "seq_region_name": "<string>",
    "start": "<integer>",
    "end": "<integer>",
    "strand": "<integer>",
    "feature_type": "<string>",
    "assembly_name": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.ensembl_get_structural_variants-parameters)

## 2) 查询 Ensembl 结构变异详情 {#dnakit.databases.ensembl_get_sv_detail}

- **作用：** 查询 Ensembl 结构变异详情。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_get_sv_detail(*, species, id=None, sv_id=None, verbose=True)`
- **输入：** 必须提供 `species`。
- **输出：** 返回字典，包含 `name`, `var_class`, `source`, `clinical_significance` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_get_sv_detail(
    species='human',
    id='esv3647486',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "name": "<string>",
  "var_class": "<string>",
  "source": "<string>",
  "clinical_significance": [
    "<string>"
  ],
  "alias": "<string>",
  "supporting_evidence": [
    {
      "name": "<string>",
      "var_class": "<string>",
      "source": "<string>",
      "clinical_significance": "<array>",
      "copy_number": "<integer>",
      "study": "<object>",
      "mappings": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.ensembl_get_sv_detail-parameters)

## 3) 按基因查询 gnomAD 结构变异 {#dnakit.databases.gnomad_get_sv_by_gene}

- **作用：** 按基因查询 gnomAD 结构变异。
- **处理方法：** 通过 `gnomADGraphQLQuery` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gnomad_get_sv_by_gene(*, gene_symbol, reference_genome=None, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** 返回字典，包含 `gene`, `dataset`, `reference_genome`, `dataset_note` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gnomad_get_sv_by_gene(
    gene_symbol='BRCA1',
    reference_genome='GRCh38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene": {
    "symbol": "<string>",
    "gene_id": "<string>",
    "chrom": "<string>",
    "start": "<integer>",
    "stop": "<integer>",
    "structural_variants": [
      "<object>"
    ]
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.gnomad_get_sv_by_gene-parameters)

## 4) 按基因组区间查询 gnomAD 结构变异 {#dnakit.databases.gnomad_get_sv_by_region}

- **作用：** 按基因组区间查询 gnomAD 结构变异。
- **处理方法：** 通过 `gnomADGraphQLQuery` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gnomad_get_sv_by_region(*, chrom, start, stop=None, reference_genome=None, end=None, verbose=True)`
- **输入：** 必须提供 `chrom`, `start`。
- **输出：** 返回字典，包含 `region`, `dataset`, `reference_genome`, `dataset_note` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gnomad_get_sv_by_region(
    chrom='17',
    start=43044295,
    stop=43070295,
    reference_genome='GRCh38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "region": {
    "chrom": "<string>",
    "start": "<integer>",
    "stop": "<integer>",
    "genes": [
      "<object>"
    ],
    "structural_variants": [
      "<object>"
    ]
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.gnomad_get_sv_by_region-parameters)

## 5) 查询 gnomAD 结构变异详情 {#dnakit.databases.gnomad_get_sv_detail}

- **作用：** 查询 gnomAD 结构变异详情。
- **处理方法：** 通过 `gnomADGraphQLQuery` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gnomad_get_sv_detail(*, variant_id=None, dataset=None, sv_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `structural_variant`, `dataset`, `reference_genome`, `dataset_note` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gnomad_get_sv_detail(
    variant_id='DEL_chr17_24e4872b',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "structural_variant": {
    "variant_id": "<string>",
    "chrom": "<string>",
    "pos": "<integer>",
    "end": "<integer>",
    "type": "<string>",
    "length": "<integer>",
    "ac": "<integer>",
    "an": "<integer>"
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_structural_variants.md#dnakit.databases.gnomad_get_sv_detail-parameters)

<!-- END GENERATED API: structural-variants -->
